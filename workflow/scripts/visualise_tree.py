from Bio import Phylo
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from os.path import join
from hashlib import sha1
import matplotlib 
import colorsys
import argparse

def string_to_color(s: str) -> str:
    if s == '?':
        return '#00DD00'
    hash_ = str(int(sha1(str(s).encode("utf-8")).hexdigest(), 16))
    color_tuple = tuple(
        (int(hash_[(i * 3):(i * 3) + 3]) % 255) / 300 for i in range(3))
    return matplotlib.colors.to_hex(color_tuple)

def darken_color(color_hex: str) -> str:
    h, l, s = matplotlib.colors.hex2color(color_hex)
    return matplotlib.colors.to_hex(colorsys.hls_to_rgb(h, .25, s))


def visualise_tree(newick_path: str, metadata_path: str, out_path: str):

    # Load metadata and styles
    tree = Phylo.read(newick_path, 'newick')
    tree.root_at_midpoint()
    tree.ladderize(reverse=True)
    sheets = pd.read_excel(metadata_path, sheet_name=None, engine='openpyxl')
    metadata = sheets['metadata'].fillna('?').astype(str).set_index('tree_name').to_dict()
    style = {}
    for sheet_name, sheet in sheets.items():
        if sheet_name == 'metadata':
            continue
        style[sheet_name] = sheet.set_index(
            sheet_name)[['color', 'marker']].to_dict(orient='index')

    show_labels = True
    n_leaves = len(tree.get_terminals())
    label_col = 'tree_new_name'

    with PdfPages(out_path) as pdf:

        for style_col in list(style) + ['region']:

            height = max(5, n_leaves / 6)
            width = 14
            fontsize = 11 - width / 5
            fig, ax = plt.subplots(1, 1, figsize=(width, height))
            ymin, ymax = (0, n_leaves)
            ax.set_ylim(ymin, ymax)

            Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)

            xmin, xmax = ax.get_xlim()

            # Get dimensions of axis in pixels
            x1, x2 = ax.get_window_extent().get_points()[:, 0]
            # Get unit scale
            xscale = (x2 - x1) / (xmax-xmin)
            # Get width of font in data units
            font_size_x_units = fontsize / xscale

            seen_style_vals = set()

            name_to_pos = {}
            name_to_style = {}

            # Update the markers and labels
            texts = [t for t in ax.texts]
            for t in texts:
                s = t.get_text().strip()

                style_val = metadata[style_col].get(s, '?')
                seen_style_vals.add(style_val)

                default_style = {'color': string_to_color(style_val), 'marker': '●'}

                if (style_val == '?') or (style_col not in style):
                    val_style = default_style
                else:
                    val_style = style[style_col].get(style_val, default_style)

                name_to_pos[s] = t.get_position()
                name_to_style[s] = val_style

                color = val_style['color']
                marker = val_style['marker']

                t.set_text(marker)
                t.set_color(color)

                t.set_size(fontsize * 1.2)
                x, y = t.get_position()
                if show_labels:
                    s = metadata[label_col].get(s, s)
                    ax.text(x + 1.8 * font_size_x_units, y, s, va='center', fontsize=fontsize)

            # Fill in the contiguous regions where the style_val is the same
            right = xmax + .15 * (xmax - xmin)
            text_right = xmax + .142 * (xmax - xmin)
            polygons = []
            to_fill = pd.DataFrame(name_to_pos).T.rename(columns={0: 'x', 1: 'y'})
            to_fill = to_fill.join(pd.DataFrame(name_to_style).T)
            to_fill['style_val'] = [metadata[style_col].get(s, '?') for s in to_fill.index]
            to_fill = to_fill.sort_values(by='y')
            poly_x = []
            poly_y: list = []
            curr_style_val = list(to_fill.style_val)[0]
            curr_color = list(to_fill.color)[0]
            for x, y, color, style_val in to_fill[['x', 'y', 'color', 'style_val']].values:
                if style_val != curr_style_val and poly_y:
                    poly_x += [right, right]
                    poly_y += [poly_y[-1], poly_y[0]]
                    polygons += [poly_x, poly_y, curr_color]
                    poly_x = []
                    poly_y = []
                    curr_style_val = style_val
                    curr_color = color
                poly_x += [x + 1.2 * font_size_x_units] * 2
                poly_y += [y - .5, y + .5]
            poly_x += [right, right]
            poly_y += [y + .5, poly_y[0]]
            polygons += [poly_x, poly_y, curr_color]
            ax.fill(*polygons, alpha=.15, ec='#555', lw=.25)

            # Add a label to each filled region
            style_val_index = 0
            curr_style_val = None
            style_val_indices = []
            for style_val in to_fill.style_val:
                if style_val != curr_style_val:
                    style_val_index += 1
                    curr_style_val = style_val
                style_val_indices.append(style_val_index)
            to_fill['style_val_index'] = style_val_indices
            mean_contiguous_y = to_fill.groupby(['style_val_index', 'style_val', 'color']).y.mean().to_frame().reset_index()
            for style_val, y, color in mean_contiguous_y[['style_val', 'y', 'color']].values:
                ax.text(text_right, y, style_val, color=darken_color(color),
                        ha='right', va='center', alpha=.5, fontsize=fontsize)

            # Make the branches thinner
            for collection in ax.collections:
                if list(collection.get_linewidths()) == [1.5]:
                    collection.set_linewidths([0.5])

            # Add a legend
            for style_val in sorted(seen_style_vals):
                default_style = {'color': string_to_color(style_val), 'marker': '●'}
                if (style_val == '?') or (style_col not in style):
                    val_style = default_style
                else:
                    val_style = style[style_col].get(style_val, default_style)
                color = val_style['color']
                scatter_style = {
                    '●': {'marker': 'o', 'fc': color, 'ec': color},
                    '■': {'marker': 's', 'fc': color, 'ec': color},
                    '○': {'marker': 'o', 'fc': '#FFF', 'ec': color},
                }[val_style['marker']]
                ax.scatter([], [], label=style_val, **scatter_style)
            ax.legend(title=' '.join(style_col.split('_')).title(), loc='upper left')

            # Hide the right and top spines
            for side in ['right', 'top']:
                ax.spines[side].set_visible(False)

            # Make the left and top axes less prominent
            faint_color = '#BBB'
            for side in ['left', 'bottom']:
                ax.spines[side].set_color(faint_color)
            ax.tick_params(axis='x', colors=faint_color)
            ax.tick_params(axis='y', colors=faint_color)
            ax.yaxis.label.set_color(faint_color)
            ax.xaxis.label.set_color(faint_color)

            ax.set_ylabel(None)

            ax.set_xlim(xmin, right)

            plt.tight_layout()
            pdf.savefig()
            plt.close()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Use metadata from a spreadsheet to color a tree visualisation')
    
    try:
        parser.add_argument('--tree', dest='newick_path', default=snakemake.input.tree)
        # convert snakemake.input.metadata to str because it comes from snakemake 
        # as a Namedlist due to the expand with organism
        parser.add_argument('--metadata', dest='metadata_path', default=str(snakemake.input.metadata))
        parser.add_argument('--out', dest='out_path', default=snakemake.output[0])
    except:
        parser.add_argument('--tree', dest='pileup_path', help='Tree in newick format', type=str, required=True)
        parser.add_argument('--metadata', dest='metadata_path', help='Metadata spreadsheet', type=str, required=True)
        parser.add_argument('--out', dest='out_path', help='Path to write pdf to', type=str, required=True)

    args = vars(parser.parse_args())
    visualise_tree(**args)
