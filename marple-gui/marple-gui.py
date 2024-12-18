import os
import cv2
import time
import shlex
import psutil
import shutil
import threading
import subprocess
import numpy as np
import tkinter as tk
from PIL import Image
from pathlib import Path
from pyzbar import pyzbar
import customtkinter as ctk
from tkinter import Menu, filedialog, ttk, messagebox

class App(ctk.CTk):
    def __init__(self):
        super().__init__()
        
        self.colmode = "light"
        self.colswitch = "#dbdbdb"
        self.update_mode()

        self.marpledir = os.path.join(os.path.join(Path.home(), 'marple'))
        self.marpleguidir = os.path.join(os.path.join(Path.home(), 'marple-gui-dev'))
        
        self.geometry("920x1020")
        self.title("MARPLE")
        self.attributes('-alpha', 0.95)
        
        self.font = ctk.CTkFont(family="Helvetica", size=16)
        self.large_font = ctk.CTkFont(family="Helvetica", size=20)
        
        self.barcode_rows = []
        self.transfer_type = 'Pgt'
        self.minknow_dir = ""
        
        # Barcode scanner variables
        self.scanner_running = False
        self.lock = threading.Lock()
        
        # Hold output lines
        self.output_lines = []
        ## Output text initialization
        self.output_text = ctk.CTkTextbox(self, width=600, height=500)
        self.output_text.configure(state="disabled")
        
        logo_path = os.path.join(self.marpleguidir, "MARPLE_logo.png")
        try:
            logo_image = ctk.CTkImage(Image.open(logo_path), size=(600, 200))
        except Exception as e:
            print(f"Error loading image: {e}")
            logo_image = None
        
        self.logo_label = ctk.CTkLabel(self, image=logo_image, text="")
        self.logo_label.pack(pady=(20, 10))

        # Menu
        self.menu_bar = Menu(self)
        self.config(menu=self.menu_bar)

        self.marple_menu = Menu(self.menu_bar, tearoff=0)
        self.menu_bar.add_cascade(label="Menu", menu=self.marple_menu)
        self.marple_menu.add_command(label="Run MARPLE", command=self.show_home)
        self.marple_menu.add_command(label="Transfer Reads", command=self.show_transfer_reads)
        self.marple_menu.add_command(label="Results", command=self.show_results)
        self.marple_menu.add_command(label="About", command=self.show_about)
        self.marple_menu.add_separator()
        self.marple_menu.add_command(label="Exit", command=self.quit)

        # Light/Dark mode toggle in the menu
        self.theme_menu = Menu(self.menu_bar, tearoff=0)
        self.menu_bar.add_cascade(label="Settings", menu=self.theme_menu)
        self.theme_menu.add_command(label=f'{"Dark mode" if self.colmode == "light" else "Light mode"}', command=self.switch_mode)

        # Run MARPLE button (Home Page)
        self.run_marple_button = ctk.CTkButton(self, text="RUN MARPLE", command=self.run_marple, corner_radius=1, font=self.large_font)
        self.run_marple_button.pack(pady=(30, 30))
        
        self.stop_button = ctk.CTkButton(self, text="STOP MARPLE", command=self.stop_marple, corner_radius=1, font=self.large_font)
        self.stop_button.pack(pady=(10, 20))
        self.snakemake_process = None

        self.transfer_in_progress = False
        self.dynamic_frame = None

        self.progress_bar = ctk.CTkProgressBar(self, orientation="horizontal")
        
        self.notification_frame = ctk.CTkFrame(self, bg_color=self.colswitch)

    def update_mode(self):
        ctk.set_appearance_mode(self.colmode)
        self.colswitch = "#dbdbdb" if self.colmode == "light" else "#2b2b2b"

    def switch_mode(self):
        self.colmode = "dark" if self.colmode == "light" else "light"
        self.update_mode()
        self.update_ui()
        # Update the theme toggle text in the menu
        self.theme_menu.entryconfig(0, label=f'{"Dark mode" if self.colmode == "light" else "Light mode"}')

    def update_ui(self):
        if self.dynamic_frame:
            self.dynamic_frame.configure(bg_color=self.colswitch)

    def show_home(self):
        self.clear_dynamic_frame()
        self.run_marple_button.pack(pady=(30, 30))
        self.stop_button.pack(pady=(10, 20))
        # self.output_text.pack(pady=(20, 20))

    def show_transfer_reads(self):
        self.clear_dynamic_frame()
        self.dynamic_frame = ctk.CTkFrame(self, bg_color=self.colswitch)
        self.dynamic_frame.pack(fill="both", expand=True, pady=20)

        # Transfer reads options (buttons already centered)
        self.select_dir_button = ctk.CTkButton(self.dynamic_frame, command=self.select_experiment, text="Select MinKNOW Directory", corner_radius=1, font=self.font)
        self.select_dir_button.pack(pady=(20, 10))

        self.expname_label = ctk.CTkLabel(self.dynamic_frame, text="Experiment Name: ", corner_radius=1, font=self.large_font)
        self.expname_label.pack(pady=(10, 20))

        self.add_row_button = ctk.CTkButton(self.dynamic_frame, text="Add Barcode Row", command=self.add_barcode_row, corner_radius=1, font=self.font)
        self.add_row_button.pack(pady=(10, 20))

        self.transfer_reads_button = ctk.CTkButton(self.dynamic_frame, command=self.transfer_reads, text="Transfer Reads", corner_radius=1, font=self.large_font)
        self.transfer_reads_button.pack(pady=(10, 20))

        # Frame to contain canvas and scrollbar
        self.canvas_frame = ctk.CTkFrame(self.dynamic_frame, bg_color=self.colswitch)
        self.canvas_frame.pack(fill="none", expand=False, padx=(5, 0), pady=5)

        self.canvas = ctk.CTkCanvas(self.canvas_frame, bg=self.colswitch, width=820, height=375)
        self.v_scrollbar = ttk.Scrollbar(self.canvas_frame, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ctk.CTkFrame(self.canvas, bg_color=self.colswitch, corner_radius=1)

        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.v_scrollbar.set)

        self.canvas.pack(side="left", fill="none", expand=True)
        self.v_scrollbar.pack(side="right", fill="y")
        
        if not self.has_barcode_rows():
            self.canvas.pack_forget()
            self.canvas_frame.pack_forget()
            self.v_scrollbar.pack_forget()

    def add_barcode_row(self):
        # Create a row container for each set of fields
        row_container = ctk.CTkFrame(self.scrollable_frame, bg_color=self.colswitch, corner_radius=1)
        row_container.pack(fill="x", anchor="center", padx=5, pady=(15,0))

        # Barcode and Sample Row
        barcode_label = ctk.CTkLabel(row_container, text="Barcode:", font=self.font)
        barcode_label.pack(side="left", padx=5)

        barcode_entry = ctk.CTkEntry(row_container, width=60, corner_radius=1, font=self.font)
        barcode_entry.pack(side="left", padx=5)
        
        # Dropdown menu for Sample Collection Barcode
        dropdown_label = ctk.CTkLabel(row_container, text="Sample Collection Barcode:", font=self.font)
        dropdown_label.pack(side="left", padx=5)

        dropdown_var = tk.StringVar(value="Recorded In-Field")
        dropdown_menu = ctk.CTkOptionMenu(row_container, corner_radius=1, variable=dropdown_var, values=["Recorded In-Field", "Extracted from Live", "No Barcode"], command=lambda choice, row=row_container: self.on_dropdown_select(choice, row))
        dropdown_menu.pack(side="left", padx=5)

        segmented_button = ctk.CTkSegmentedButton(row_container, values=['Pgt', 'Pst'], command=lambda choice, row=row_container: self.on_segmented_button_click(choice, row), corner_radius=1, font=self.font)
        segmented_button.set('Pgt')
        segmented_button.pack(side="left", padx=(5, 10))

        remove_button = ctk.CTkButton(row_container, text="Remove", command=lambda: self.remove_barcode_row(row_container), corner_radius=1, font=self.font)
        remove_button.pack(side="left", padx=5)

        # Metadata container for the dropdown selection
        metadata_container = ctk.CTkFrame(self.scrollable_frame, bg_color=self.colswitch, corner_radius=1)
        metadata_container.pack(fill="x", anchor="center", padx=5)

        # Add row and metadata fields to barcode_rows for tracking
        self.barcode_rows.append({
            "row_container": row_container,
            "barcode": barcode_entry,
            "transfer_type": segmented_button,
            "dropdown_var": dropdown_var,
            "dropdown_menu": dropdown_menu,
            "metadata_container": metadata_container
        })

        # Show default option's fields by default
        self.on_dropdown_select("Recorded In-Field", row_container)

        # Update scrollable region
        self.scrollable_frame.update_idletasks()
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))

        # Show scrollbar only if rows are present
        if len(self.barcode_rows) > 0:
            self.v_scrollbar.pack(side="right", fill="y")
            self.canvas.pack(fill="none", expand=False, padx=(5, 0), pady=20)
            self.canvas_frame.pack(fill="none", expand=False, padx=(5, 0), pady=5)

    def on_dropdown_select(self, choice, row):
        # Find the metadata container for the row
        metadata_container = next(item["metadata_container"] for item in self.barcode_rows if item["row_container"] == row)
        self.clear_metadata_fields(metadata_container)
        for item in self.barcode_rows:
            if item["row_container"] == row and "metadata_container2" in item:
                item["metadata_container2"].destroy()
                del item["metadata_container2"]
                break

        if choice == "Recorded In-Field":
            self.scan_barcode_option(metadata_container, row)
        elif choice == "Extracted from Live":
            self.scan_two_barcodes_option(metadata_container, row)
        elif choice == "No Barcode":
            self.show_no_barcode_option(metadata_container, row)
            
    def scan_barcode_cam(self, marple_barcode=False):
        cam_index = self.find_external_camera()
        if cam_index is None:
            messagebox.showerror("Error", "No camera found!")
            return
        
        cap = cv2.VideoCapture(cam_index)
        if not cap.isOpened():
            messagebox.showerror("Error", "Could not open any camera.")
            return None
        
        barcode_data = None
        self.scanner_running = True

        # Add timeout to limit cpu/gpu usage
        start_time = time.time()

        while self.scanner_running:
            elapsed_time = time.time() - start_time
            if elapsed_time > 20:  # 20-second timeout
                messagebox.showwarning("Warning","Camera timeout. No barcode detected within 20 seconds.")
                break

            ret, frame = cap.read()
            if not ret:
                self.printin("Error: Failed to capture image.")
                break

            # Trying to increase the chances of successful scan by enhancing image
            # Convert frame to grayscale
            gray_frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

            # Apply a sharpening filter
            kernel = np.array([[0, -1, 0], 
                            [-1, 5, -1], 
                            [0, -1, 0]])
            
            # Combine sharpening with contrast enhancement using CLAHE (Contrast Limited Adaptive Histogram Equalization))
            
            clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
            equalized = clahe.apply(gray_frame)

            # Sharpen after histogram equalization
            sharp_frame = cv2.filter2D(equalized, -1, kernel)

            # Adjust brightness and contrast (clip values to stay in valid range)
            alpha = 2    # Contrast control (1.0-3.0)
            beta = 20    # Brightness control (0-100)
            enhanced_frame = cv2.convertScaleAbs(sharp_frame, alpha=alpha, beta=beta)

            # Upscale image for better decoding
            upscale_factor = 2
            enhanced_frame = cv2.resize(enhanced_frame, 
                                        (frame.shape[1] * upscale_factor, frame.shape[0] * upscale_factor), 
                                        interpolation=cv2.INTER_LINEAR)

            barcodes = pyzbar.decode(enhanced_frame)

            for barcode in barcodes:
                barcode_data = barcode.data.decode("utf-8")
                if marple_barcode and barcode_data.startswith("M"):
                    self.scanner_running = False
                    messagebox.showinfo("Info",f"MARPLE Barcode detected: {barcode_data}")
                    break
                elif marple_barcode and not barcode_data.startswith("M"):
                    # I used to be able to just continue scanning until a valid barcode was detected
                    # but now it seems to be stuck in an infinite loop if an invalid barcode is scanned
                    # so I'm breaking out of the loop and returning None
                    self.scanner_running = False
                    messagebox.showerror("Error","Invalid barcode. Please scan a MARPLE barcode.")
                    break
                else:
                    self.printin(f"Barcode detected: {barcode_data}")
                    self.scanner_running = False
                    break

            cv2.imshow("Scan Barcode", frame)

            if (cv2.waitKey(1) & 0xFF == ord('q')) or (cv2.getWindowProperty("Scan Barcode", cv2.WND_PROP_VISIBLE) < 1):
                self.scanner_running = False
                break

        cap.release()
        cv2.destroyAllWindows()
        return barcode_data

    def find_external_camera(self):
        max_inputs = 8
        camera_indices = []
        
        for i in range(max_inputs):
            cap = cv2.VideoCapture(i)
            if cap.isOpened():
                # Test if camera works
                ret, _ = cap.read()
                if ret:
                    camera_indices.append(i)
                cap.release()
        
        # Identify external webcam (assume external by resolution)
        for index in camera_indices:
            cap = cv2.VideoCapture(index)
            if cap.isOpened():
                # Retrieve resolution
                width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
                height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
                print(f"Camera {index}: {width}x{height}")
                cap.release()
                if width > 1280 or height > 720:
                    return index
        
        if camera_indices[0] == 0:
            self.printin(f"External webcam not found. Falling back to the default camera ({height}p).")

        # Default to the first available camera if no external is found
        return camera_indices[0] if camera_indices else None

    def toggle_scanner(self, container, row, marple_toggle=False):
        self.run_scanner(container, row, marple_toggle)

    def run_scanner(self, container, row, marple_toggle=False):
        with self.lock:
            self.scanner_running = True
        if marple_toggle:
            marple_barcode_data = self.scan_barcode_cam(marple_barcode=True)
            self.after(0, self.update_marple_entry, container, row, marple_barcode_data)
        else:
            odk_barcode_data = self.scan_barcode_cam(marple_barcode=False)
            self.after(0, self.update_odk_entry, container, row, odk_barcode_data)
        with self.lock:
            self.scanner_running = False

    def update_marple_entry(self, container, row, marple_barcode_data):
        with self.lock:
            for item in self.barcode_rows:
                if item["row_container"] == row:
                    item.update({
                        "marple_barcode": marple_barcode_data
                    })
                    break

    def update_odk_entry(self, container, row, odk_barcode_data):
        with self.lock:
            for item in self.barcode_rows:
                if item["row_container"] == row:
                    item.update({
                        "odk_barcode": odk_barcode_data
                    })
                    break

    def scan_barcode_option(self, container, row):
        barcode_label = ctk.CTkButton(container, text="Scan MARPLE Barcode", font=self.font, command=lambda: self.toggle_scanner(container, row, marple_toggle=True), corner_radius=1)
        barcode_label.pack(side="left", padx=5)

    def scan_two_barcodes_option(self, container, row):
        qrcode1_label = ctk.CTkButton(container, text="Scan ODK Barcode", font=self.font, command=lambda: self.toggle_scanner(container, row, marple_toggle=False), corner_radius=1)
        qrcode1_label.pack(side="left", padx=5)

        qrcode2_label = ctk.CTkButton(container, text="Scan MARPLE Barcode", font=self.font, command=lambda: self.toggle_scanner(container, row, marple_toggle=True), corner_radius=1)
        qrcode2_label.pack(side="left", padx=5)

    def show_no_barcode_option(self, container, row):
        # Add sample metadata entries
        sample = ctk.CTkEntry(container, width=200, placeholder_text="Sample Name", corner_radius=1, font=self.font)
        sample.pack(side="left", padx=2)

        collection_date = ctk.CTkEntry(container, width=200, placeholder_text="Collection Date", font=self.font, corner_radius=1)
        collection_date.pack(side="left", padx=2)

        location = ctk.CTkEntry(container, width=200, placeholder_text="Location", font=self.font, corner_radius=1)
        location.pack(side="left", padx=2)

        country = ctk.CTkEntry(container, width=200, placeholder_text="Country", font=self.font, corner_radius=1)
        country.pack(side="left", padx=2)
        
        metadata_container2 = ctk.CTkFrame(self.scrollable_frame, bg_color=self.colswitch, corner_radius=1)
        metadata_container2.pack(fill="x", anchor="center", padx=5)
        
        collector_name = ctk.CTkEntry(metadata_container2, width=200, placeholder_text="Collector's Name", font=self.font, corner_radius=1)
        collector_name.pack(side="left", padx=2)
        
        cultivar = ctk.CTkEntry(metadata_container2, width=200, placeholder_text="Cultivar", font=self.font, corner_radius=1)
        cultivar.pack(side="left", padx=2)

        treat = ctk.CTkEntry(metadata_container2, width=200, placeholder_text="Treatment", font=self.font, corner_radius=1)
        treat.pack(side="left", padx=2)
        
        # Update the row's metadata fields
        for item in self.barcode_rows:
            if item["metadata_container"] == container:
                item.update({
                    "sample": sample,
                    "collection_date": collection_date,
                    "location": location,
                    "country": country,
                    "metadata_container2": metadata_container2,
                    "collector_name": collector_name,
                    "cultivar": cultivar,
                    "treat": treat
                })
                break

    def clear_metadata_fields(self, container):
        for widget in container.winfo_children():
            widget.pack_forget()            

    def remove_barcode_row(self, row):
        # Find the index of the row to remove
        index_to_remove = None
        for index, barcode_row in enumerate(self.barcode_rows):
            if barcode_row["row_container"] == row:
                index_to_remove = index
                break

        if index_to_remove is not None:
            # Remove the row from the barcode_rows list
            del self.barcode_rows[index_to_remove]

            # Destroy the row widget and metadata container
            row.destroy()
            metadata_container = barcode_row["metadata_container"]
            self.clear_metadata_fields(metadata_container)
            metadata_container.destroy()

            # Remove metadata_container2 if it exists
            if "metadata_container2" in barcode_row:
                barcode_row["metadata_container2"].destroy()

            # Update the canvas and scrollbar if necessary
            self.scrollable_frame.update_idletasks()
            self.canvas.configure(scrollregion=self.canvas.bbox("all"))

            if not self.has_barcode_rows():
                self.canvas.pack_forget()
                self.canvas_frame.pack_forget()
                self.v_scrollbar.pack_forget()
    
    def has_barcode_rows(self):
        return len(self.barcode_rows) > 0

    def on_segmented_button_click(self, choice, row):
        # Find the index of the row to update
        index = next((index for index, item in enumerate(self.barcode_rows) if item["row_container"] == row), None)
        if index is not None:
            # Update the transfer type for that specific row
            self.barcode_rows[index]["transfer_type"].set(choice)

    def select_experiment(self):
        self.minknow_default_dir = "/var/lib/minknow/data"
        self.minknow_dir = filedialog.askdirectory(initialdir=self.minknow_default_dir)
        if self.minknow_dir and (self.minknow_dir != self.minknow_default_dir):
            expname = os.path.basename(self.minknow_dir)
            self.expname_label.configure(text=f"Experiment Name: {expname}")

    def transfer_reads(self):
        if not self.minknow_dir or self.minknow_dir == self.minknow_default_dir:
            messagebox.showerror("Error","MinKNOW directory not selected.")
            return

        for row in self.barcode_rows:
            # Check for each field depending on dropdown selection
            if row["dropdown_var"].get() == "Recorded In-Field":
                marple_barcode = row.get("marple_barcode")
                if not marple_barcode or not marple_barcode.strip().startswith("M"):
                    messagebox.showerror("Please scan a valid QR code.")
                    return

            elif row["dropdown_var"].get() == "Extracted from Live":
                odk_barcode = row.get("odk_barcode")
                marple_barcode = row.get("marple_barcode")
                if not odk_barcode or not marple_barcode or not marple_barcode.strip().startswith("M"):
                    messagebox.showerror("Please scan valid QR codes.")
                    return

            elif row["dropdown_var"].get() == "No Barcode":
                if any(not field.get().strip() for field in [
                    row["sample"], row["collection_date"], row["location"],
                    row["country"], row["collector_name"], row["cultivar"]]):
                    messagebox.showerror("Please fill in all metadata fields.")
                    
        # Check if a transfer is already in progress
        if self.transfer_in_progress:
            messagebox.showwarning("Warning","Transfer Reads Process Running. Wait until it finishes.")
            return

        # Start the transfer process
        self.transfer_in_progress = True
        self.progress_bar.pack(pady=(10, 20))
        self.progress_bar.configure(mode="indeterminate")
        self.progress_bar.start()

        # Run the task in a separate thread to avoid freezing the UI
        self.transfer_thread = threading.Thread(target=self.process_reads)
        self.transfer_thread.start()
            
    def process_reads(self):
        try:
            total_barcodes = len(self.barcode_rows)
            if total_barcodes == 0:
                messagebox.showerror("Error","No barcodes to process.")
                return

            # Read existing metadata
            metadata_file = os.path.join(self.marpledir, 'sample_metadata.csv')
            if not os.path.exists(metadata_file):
                with open(metadata_file, 'w') as f:
                    f.write('Experiment,Barcode,SampleName,Pathogen,CollectionDate,CollectorsName,Location,Country,Cultivar,Treatment,ODKcode\n')

            existing_metadata = []
            if os.path.exists(metadata_file):
                with open(metadata_file, 'r') as f:
                    existing_metadata = f.readlines()

            for row in self.barcode_rows:
                experiment = os.path.basename(self.minknow_dir)
                barcode_entry = row["barcode"]
                sample_entry = row.get("sample") or row.get("marple_barcode")
                segmented_button = row["transfer_type"]
                meta_date = row.get("collection_date")
                meta_name = row.get("collector_name")
                meta_loc = row.get("location")
                meta_country = row.get("country")
                meta_cultivar = row.get("cultivar")
                meta_treat = row.get("treat")
                meta_odk = row.get("odk_barcode")

                barcode = format(int(barcode_entry.get().strip()), '02d')
                if hasattr(sample_entry, "strip"):
                    sample = sample_entry.strip() if sample_entry else ""
                else:
                    sample = sample_entry.get() if sample_entry else ""
                pathogen = segmented_button.get()  # Get the transfer type from the segmented button

                if not barcode:
                    continue

                # Check if the sample and pathogen already exist in the metadata
                new_metadata_line = f"{experiment},{barcode},{sample},{pathogen},{meta_date.get().strip() if meta_date else ''},{meta_name.get().strip() if meta_name else ''},{meta_loc.get().strip() if meta_loc else ''},{meta_country.get().strip() if meta_country else ''},{meta_cultivar.get().strip() if meta_cultivar else ''},{meta_treat.get().strip() if meta_treat else ''},{meta_odk.strip() if meta_odk else ''}\n"
                existing_metadata = [line for line in existing_metadata if not (sample in line and pathogen in line)]

                for root, dirs, files in os.walk(self.minknow_dir):
                    for dir in dirs:
                        if dir == 'pass':
                            barcode_dir = os.path.join(root, dir, f'barcode{barcode}')
                            if os.path.exists(barcode_dir):
                                try:
                                    output_file = os.path.join(self.marpledir, 'reads', pathogen.lower(), f'{sample}.fastq.gz')
                                    command = f"cat {barcode_dir}/*.fastq.gz > {output_file}"
                                    subprocess.run(command, shell=True)
                                    self.printin(f"Successfully transferred reads for barcode{barcode} to {output_file}")

                                    existing_metadata.append(new_metadata_line)

                                except Exception as e:
                                    messagebox.showerror("Error",f"An error occurred while processing barcode{barcode}: {e}")
                            else:
                                messagebox.showerror("Error",f"Barcode{barcode} not found in MinKNOW directory.")

            # Write the updated metadata back to the file
            with open(metadata_file, 'w') as f:
                f.writelines(existing_metadata)

            backup_dir = '/marple/upload'
            if not os.path.exists(backup_dir):
                os.makedirs(backup_dir)
            shutil.copy(metadata_file, os.path.join(backup_dir, 'sample_metadata.csv'))

        finally:
            self.transfer_in_progress = False
            self.stop_progress_bar()
    
    def show_about(self):
        self.clear_dynamic_frame()
        self.dynamic_frame = ctk.CTkFrame(self, bg_color=self.colswitch)
        self.dynamic_frame.pack(fill="both", expand=True, pady=20)

        about_label = ctk.CTkLabel(self.dynamic_frame, text="MARPLE Diagnostics:\npoint-of-care, strain-level disease diagnostics and\nsurveillance tool for complex fungal pathogens\n\nVersion: 2.0-alpha", font=self.large_font)
        about_label.pack(pady=(20, 20))

        # devs_label = ctk.CTkLabel(self.dynamic_frame, text="Software Development:\nLoizos Savva")
        # devs_label.pack(pady=(20, 20))
        
        copyright_label = ctk.CTkLabel(self.dynamic_frame, text="© 2024 Saunders Lab")
        copyright_label.pack()

    def stop_progress_bar(self):
        self.progress_bar.stop()
        self.progress_bar.pack_forget()
        
    def clear_dynamic_frame(self):
        if self.dynamic_frame:
            self.dynamic_frame.destroy()
        
        self.run_marple_button.pack_forget()
        self.stop_button.pack_forget()
        self.output_text.pack_forget()
            
    def run_marple(self):
        # Start progress bar
        self.progress_bar.pack(pady=(10, 20))
        self.progress_bar.configure(mode="indeterminate")
        self.progress_bar.start()

        # Ensure output_text exists
        if not hasattr(self, 'output_text'):
            self.output_text = ctk.CTkTextbox(self, width=600, height=400)
            self.output_text.pack(pady=(20, 20))
            self.output_text.configure(state="disabled")

        # Run the Snakemake process in a separate thread
        thread = threading.Thread(target=self.run_snakemake)
        thread.start()
        # thread.join()
        
        self.output_text.pack(pady=(20, 20))

    def run_snakemake(self):
        if shutil.which("mamba") is not None:
            try:
                env_list = subprocess.check_output(["mamba", "env", "list"], text=True)
                if "marple-env" in env_list:
                    self.stop_marple(forced=False)
                    self.start_marple()
                else:
                    messagebox.showerror("Error","mamba environment marple-env not found.")
            except subprocess.CalledProcessError as e:
                messagebox.showerror("Error",f"Error running command: {e}")
        else:
            messagebox.showerror("Error","mamba not found.")
        
        self.progress_bar.stop()
        self.progress_bar.pack_forget()

    def start_marple(self):
        try:
            # Start the Snakemake process with unbuffered output
            self.snakemake_process = subprocess.Popen(
                ["snakemake", "--cores", "all", "--rerun-incomplete"],
                cwd=self.marpledir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, bufsize=1, universal_newlines=True
            )
            
            # Use a thread to capture output and display it in real-time
            threading.Thread(target=self.capture_output, args=(self.snakemake_process.stdout,), daemon=True).start()
            threading.Thread(target=self.capture_output, args=(self.snakemake_process.stderr,), daemon=True).start()
            
            self.printin("MARPLE Snakemake workflow started.")
        except Exception as e:
            messagebox.showerror("Error",f"Failed to start Snakemake: {e}")

    def capture_output(self, stream):
        # Append captured output to the output_lines
        for line in iter(stream.readline, ''):
            self.output_lines.append(line)  # Store captured output
            self.output_text.configure(state="normal")
            self.output_text.insert("end", line)
            self.output_text.see("end")
            self.output_text.configure(state="disabled")
        stream.close()

    def stop_marple(self, forced=True):
        # Stop the Snakemake process if running
        if self.snakemake_process and self.snakemake_process.poll() is None:
            self.snakemake_process.terminate()
            messagebox.showinfo("MARPLE process terminated.")
            
            # Kill child processes of Snakemake (if any)
            try:
                parent = psutil.Process(self.snakemake_process.pid)
                children = parent.children(recursive=True)  # Get child processes
                for child in children:
                    self.printin(f"Terminating child process: {child.name()} (PID {child.pid})")
                    child.terminate()
                    child.wait()
                self.printin("All child processes terminated.")
            except Exception as e:
                self.printin(f"Failed to terminate child processes: {e}")
        
        # Check for independent processes like FastTree
        try:
            for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
                cmdline = proc.info['cmdline']
                if cmdline and 'fasttree' in cmdline:
                    self.printin(f"Terminating FastTree process (PID {proc.info['pid']})")
                    proc.terminate()
                    proc.wait()
        except Exception as e:
            self.printin(f"Error terminating FastTree: {e}")
        
        self.unlock_snakemake()
        
        if forced:
            self.output_text.pack_forget()

    def unlock_snakemake(self):
        try:
            subprocess.run(["mamba", "run", "-n", "marple-env", "snakemake", "--unlock"], cwd=self.marpledir)
        except Exception as e:
            self.printin(f"Failed to unlock Snakemake: {e}")

    def toggle_output_display(self):
        if self.output_text.winfo_viewable():
            self.output_text.pack_forget()
        else:
            self.output_text.pack(pady=(20, 20))
    
    def check_process(self, process, log_file):
        if process.poll() is None:
            self.after(100, self.check_process, process, log_file)
        else:
            self.progress_bar.stop()
            self.progress_bar.pack_forget()
            status = process.returncode
            if status != 0:
                self.printin(f"Command failed with exit status {status}. Error log:")
                with open(log_file, "r") as log:
                    self.printin(log.read())
            else:
                self.printin("MARPLE run completed successfully.")

    def printin(self, stdout):
        self.notification_frame.pack(side="top", fill="both", expand="True")
        for widget in self.notification_frame.winfo_children():
            widget.destroy()
        
        label = ctk.CTkLabel(self.notification_frame, text=stdout, bg_color=self.colswitch)
        label.pack(side="top", fill="both", expand="True")
        
        self.after(5000, self.clear_notification)
    
    def clear_notification(self):
        for widget in self.notification_frame.winfo_children():
            widget.destroy()
        self.notification_frame.pack_forget()

    def show_results(self):
        self.clear_dynamic_frame()
        self.dynamic_frame = ctk.CTkFrame(self, bg_color=self.colswitch)
        self.dynamic_frame.pack(fill="both", expand=True, pady=20)

        # Two columns
        self.dynamic_frame.grid_columnconfigure(0, weight=1)
        self.dynamic_frame.grid_columnconfigure(1, weight=1)

        # Pgt Column
        pgt_label = ctk.CTkLabel(self.dynamic_frame, text="Pgt", font=self.large_font, anchor="center")
        pgt_label.grid(row=0, column=0, pady=(10, 10))

        pgt_tree_button = ctk.CTkButton(self.dynamic_frame, text="Tree", command=lambda: self.open_file("pgt", "trees", "pgt_all.pdf"), corner_radius=1, font=self.font)
        pgt_tree_button.grid(row=1, column=0, pady=(10, 10))

        pgt_report_button = ctk.CTkButton(self.dynamic_frame, text="Report", command=lambda: self.open_file("pgt", "report", "pgt.multiqc.html"), corner_radius=1, font=self.font)
        pgt_report_button.grid(row=2, column=0, pady=(10, 10))

        # Pst Column
        pst_label = ctk.CTkLabel(self.dynamic_frame, text="Pst", font=self.large_font, anchor="center")
        pst_label.grid(row=0, column=1, pady=(10, 10))

        pst_tree_button = ctk.CTkButton(self.dynamic_frame, text="Tree", command=lambda: self.open_file("pst", "trees", "pst_all.pdf"), corner_radius=1, font=self.font)
        pst_tree_button.grid(row=1, column=1, pady=(10, 10))

        pst_report_button = ctk.CTkButton(self.dynamic_frame, text="Report", command=lambda: self.open_file("pst", "report", "pst.multiqc.html"), corner_radius=1, font=self.font)
        pst_report_button.grid(row=2, column=1, pady=(10, 10))

    def open_file(self, type, subdir, filename):
        file_path = os.path.join(self.marpledir, "results", type, subdir, filename)
        if filename.endswith(".pdf"):
            print(f'xdg-open {file_path}')
            os.system(f"xdg-open {file_path}")
        elif filename.endswith(".html"):
            os.system(f"xdg-open {file_path}")
        

app = App()
app.mainloop()
