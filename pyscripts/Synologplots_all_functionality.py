#!/bin/python

import argparse, cairo, math, os, sys
from collections import defaultdict
import numpy as np


def get_arguments():
    """The name is pretty self explantory"""

    arg_description = (
        """ Trying to recreate a synolog plot. Hopefully it all goes well """
    )

    parser = argparse.ArgumentParser(description=arg_description)

    parser.add_argument(
        "-sA", "--speciesA", help="Species One Name", default=None, type=str
    )

    parser.add_argument(
        "-sB", "--speciesB", help="Species Two Name", default=None, type=str
    )

    parser.add_argument(
        "-sC", "--speciesC", help="Species Three Name", default=None, type=str
    )

    parser.add_argument(
        "-S",
        "--synolog_dir",
        help="Path to Synolog Output",
        required=True,
        type=str,
    )

    parser.add_argument(
        "-c",
        "--color_scheme",
        help="Color Scheme: 1 (Chromosome), 2 (Synteny Cluster)",
        default=1,
        type=int,
    )

    args = parser.parse_args()
    spA = args.speciesA
    spB = args.speciesB
    spC = args.speciesC
    synolog_dir = args.synolog_dir
    color_scheme = args.color_scheme

    return spA, spB, spC, synolog_dir, color_scheme


###############
### CLASSES ###
###############

# loop through and make each row an object to store in a list
class Single_Ortho_Row:
    def __init__(self, spA, spA_chrom, spB, spB_chrom, startA, endA, startB, endB):
        self.spA = spA
        self.spB = spB
        self.spA_chrom = spA_chrom
        self.spB_chrom = spB_chrom
        self.start_posA = startA
        self.end_posA = endA
        self.start_posB = startB
        self.end_posB = endB


# store color params for cairo.context.set_source_rgba()
class Chrom_Colors:
    def __init__(self, red, green, blue, alpha):
        self.red = red
        self.green = green
        self.blue = blue
        self.alpha = alpha


# this was a brilliant idea from Angel, so I am going to use it as well
class Pairing_Colors:
    def __init__(self, value):
        colors_dict = {
            0: [0.75, 0.25, 0.25, 0.8],  # light red
            1: [0.7, 0.2, 0.7, 0.8],  # pink
            2: [0.9, 0.9, 0.2, 0.8],  # yellow
            3: [0.2, 0.8, 0.2, 0.8],  # lime
            4: [0.1, 0.1, 0.9, 0.8],  # blue
            5: [0.2, 0.5, 0.2, 0.8],  # green
            6: [0.6, 0.2, 0.9, 0.8],  # purple
            7: [0.2, 0.5, 0.5, 0.8],  # light navy-ish
        }

        colors = colors_dict[value % len(colors_dict)]
        self.red = colors[0]
        self.blue = colors[1]
        self.green = colors[2]
        self.alpha = colors[3]


#########################
### Parsing Functions ###
#########################


def get_scaffold_lengths(synolog_dir, spA, spB):
    """parse the chromosomes file for chrom lengths"""

    path2chrom = os.path.join(synolog_dir, "chromosomes.tsv")

    spA_lengths = {}
    spB_lengths = {}

    with open(path2chrom, "r") as chrom:
        for line in chrom:
            chrom_fields = line.strip().split("\t")
            # only keep scaffolds atleast 1mb in length
            if chrom_fields[1] == spA:
                if int(chrom_fields[3]) >= 1e6:
                    spA_lengths[chrom_fields[2]] = int(chrom_fields[3])
            elif chrom_fields[1] == spB:
                if int(chrom_fields[3]) >= 1e6:
                    spB_lengths[chrom_fields[2]] = int(chrom_fields[3])

    # sort to get scaffolds longest to shortest
    spA_dict = {
        scaf: scaf_len
        for scaf, scaf_len in sorted(
            spA_lengths.items(), key=lambda item: item[1], reverse=True
        )
    }

    spB_dict = {
        scaf: scaf_len
        for scaf, scaf_len in sorted(
            spB_lengths.items(), key=lambda item: item[1], reverse=True
        )
    }

    return spA_dict, spB_dict


def get_ortho_rows(spA, spB, synolog_dir):
    """parse the membership files to get the relevant information per row"""

    membership_file = f"{spA}-{spB}_consyn_cluster_membership.tsv"
    path2member = os.path.join(synolog_dir, membership_file)

    # create a list to store all the Single_Ortho_Row objects
    ortho_matches = []

    with open(path2member, "r") as syno:
        for line in syno:
            if line.startswith("#"):
                continue
            syno_fields = line.strip().split("\t")
            spA_name = syno_fields[2]
            spB_name = syno_fields[3]
            spA_chrom = syno_fields[6]
            spB_chrom = syno_fields[12]
            spA_start = int(syno_fields[8])
            spA_end = int(syno_fields[9])
            spB_start = int(syno_fields[14])
            spB_end = int(syno_fields[15])
            ortho_row = Single_Ortho_Row(
                spA_name,
                spA_chrom,
                spB_name,
                spB_chrom,
                spA_start,
                spA_end,
                spB_start,
                spB_end,
            )
            ortho_matches.append(ortho_row)

    return ortho_matches


############################
### ORGANIZING FUNCTIONS ###
############################


def quick_color_generator(current_iteration):
    """will be used to alternate colors based on if iteration is even/odd"""

    iter_remainder = current_iteration % 2
    if iter_remainder == 0:
        chrom_color = Chrom_Colors(0.576, 0.914, 0.745, 1)  # sea foam
    else:
        chrom_color = Chrom_Colors(0.7, 0.7, 0.7, 1)  # light grey

    return chrom_color


def normalize_positions(starting_x, x1, x2, scaffold_length, genome_size, image_width):
    """Attempting to scale down the coordinates to pdf dimensions"""

    x = int((x1 + x2) / 2)  # midpoint

    # get relative length of scaffold to full genome
    relative_scaffold_len = scaffold_length / genome_size

    # find what proportion of the allocation image width space that length takes
    occupied_length = relative_scaffold_len * image_width

    # see where the midpoint lands in that chrom/scaffold
    midpoint_location = x / scaffold_length

    # use the starting location to find its scaled location
    new_pos = starting_x + (midpoint_location * occupied_length)

    return math.ceil(new_pos)


def merge_dicts(spA_dict, spB_dict, ortho_rows):
    """function to put all ortho rows into a dict where keys are scaffold matches"""
    full_dict = defaultdict(list)

    for ortho_row in ortho_rows:
        spA_chrom = ortho_row.spA_chrom
        spB_chrom = ortho_row.spB_chrom
        # check if this comparison survived 1mb filter
        if spA_chrom in spA_dict.keys():
            if spB_chrom in spB_dict.keys():
                spA_chrom_length = spA_dict[spA_chrom]
                spB_chrom_length = spB_dict[spB_chrom]
                spA_start = ortho_row.start_posA
                spA_end = ortho_row.end_posA
                spB_start = ortho_row.start_posB
                spB_end = ortho_row.end_posB
                full_dict[(spA_chrom, spB_chrom)].append(
                    [
                        spA_chrom_length,
                        spB_chrom_length,
                        spA_start,
                        spA_end,
                        spB_start,
                        spB_end,
                    ]
                )

    return full_dict


def chrom_reorder(sorted_spA_dict, sorted_spB_dict, ortho_rows):
    """Sort the scaffolds/chroms on the bottom based on top scaffold/chrom order"""

    # these scaffolds will be plotted first to align with spA
    first_chromsB_lists = {}

    for spA_chrom in sorted_spA_dict.keys():
        spB_chrom_ranges = defaultdict(int)
        for ortho_row in ortho_rows:
            if spA_chrom == ortho_row.spA_chrom:
                spB_chrom = ortho_row.spB_chrom
                # since some chroms are filtered out (< 1mb), need to do a check
                if spB_chrom in sorted_spB_dict.keys():
                    spB_start = ortho_row.start_posB
                    spB_end = ortho_row.end_posB
                    spB_range = spB_end - spB_start  # length of covered region
                    spB_chrom_ranges[spB_chrom] += spB_range
        # filter to the most longest region
        longest_region = 0
        longest_spB_chrom = None
        for spB_chrom, spB_region in spB_chrom_ranges.items():
            if spB_region > longest_region:
                longest_region = spB_region
                longest_spB_chrom = spB_chrom
        if longest_spB_chrom is None:
            continue
        else:
            first_chromsB_lists[longest_spB_chrom] = sorted_spB_dict[longest_spB_chrom]

    return first_chromsB_lists


def establish_starting_positions(sp_dict, genome_size, image_width, starting_pos):
    """use a sorted dictionary to determine the starting position for each scaffold"""

    # size of allocated x-axis
    allocated_sites = image_width / genome_size

    scaled_sp_dict = {}

    for scaffold_name, scaffold_size in sp_dict.items():
        # assign 1st, then see how much is occupied for the next items
        scaled_sp_dict[scaffold_name] = starting_pos
        scaled_size = allocated_sites * scaffold_size
        occupied_size = starting_pos + scaled_size
        starting_pos = math.ceil(occupied_size + 1)

    return scaled_sp_dict


def check_for_reverse(ortho_rows):
    """Apply a linear agression to determine if scaffolds should be reversed"""
    reverse = None

    # get possitions
    spA_list = []
    spB_list = []

    # add start positions
    for rows in ortho_rows:
        spA_start = rows[2]
        spB_start = rows[4]
        spA_list.append(spA_start)
        spB_list.append(spB_start)

    # get average
    spA_avg = np.mean(spA_list)
    spB_avg = np.mean(spB_list)

    avg_multiplied = 0
    avg_squared = 0

    # num of positions
    pos_ct = len(ortho_rows)

    for i in range(pos_ct):
        posA = spA_list[i]  # get values
        posB = spB_list[i]
        avg_multiplied += posA * posB  # multiply & add
        avg_squared += posA ** 2  # square for spA

    # find averages
    avg_multiplied = avg_multiplied / pos_ct
    avg_squared = avg_squared / pos_ct

    # get slope
    slope = ((spA_avg * spB_avg) - avg_multiplied) / ((spA_avg ** 2) - avg_squared)

    if slope >= 0:
        reverse = False
    else:
        reverse = True

    return reverse


def Reverse_Ortho_Rows(ortho_rows):
    """Reverse the positions to plot in reverse"""

    new_ortho_rows = []

    # get num of rows to get opposite index
    total_rows = len(ortho_rows) - 1

    for i in range(total_rows):
        reversed_i = total_rows - i
        row_i = ortho_rows[i]
        row_r = ortho_rows[reversed_i]  # opposite position
        spA_genome_length = row_i[0]
        spB_genome_length = row_r[1]
        spA_start = row_i[2]
        spA_end = row_i[3]
        spB_start = row_r[4]
        spB_end = row_r[5]
        new_ortho_rows.append(
            [
                spA_genome_length,
                spB_genome_length,
                spA_start,
                spA_end,
                spB_start,
                spB_end,
            ]
        )

    return new_ortho_rows


def create_ticks(genome_size):
    """Find the appropriate marker lengths for tick marks"""

    marker_sets = {
        1e3: "K",
        1e4: "K",
        1e5: "K",  # thousands
        1e6: "M",
        1e7: "M",
        1e8: "M",  # millions
        1e9: "B",
        1e10: "B",
        1e11: "B",  # billions
    }

    for marker_length in marker_sets.keys():
        if genome_size < marker_length:
            upper_limit = list(marker_sets.keys()).index(marker_length)
            # get the previous marker (i.e., size down)
            marker = list(marker_sets.values())[upper_limit - 1]
            break

    # check how many ticks can fit given a genome size
    # multiplier will be used to show increments in bp length
    if marker == "K":
        tick_ct = math.floor(genome_size / 5e4)
        multiplier = 50
    elif marker == "M":
        tick_ct = math.floor(genome_size / 5e7)
        multiplier = 50
    elif marker == "B":
        marker = "M"
        tick_ct = math.floor(genome_size / 5e7)
        multiplier = 500
    return marker, tick_ct, multiplier


##########################
### PLOTTING FUNCTIONS ###
##########################


def Draw_Ticks(ctx, genome_length, genome_orientation, marker_direction, image_width):
    """Draw the tick markers given a genome's size and top/bottom orientation"""

    # need to offset from the chrom lines
    # genome_orientation == y pos of top/bottom chrom
    if marker_direction == "top":
        mark_ypos1 = genome_orientation - 23
        mark_ypos2 = genome_orientation - 75
        label_ypos = genome_orientation - 80
    elif marker_direction == "bottom":
        mark_ypos1 = genome_orientation + 23
        mark_ypos2 = genome_orientation + 75
        label_ypos = genome_orientation + 115

    marker, tick_ct, multiplier = create_ticks(genome_length)
    allocated_pos = image_width / tick_ct

    for i in range(1, tick_ct + 1):
        marker_num = i * multiplier  # interval scheme
        marker_label = f"{str(marker_num)}{marker}"
        marker_pos = 100 + (allocated_pos * i)  # 100 pos offset for image

        # plot each tick and its label
        ctx.move_to(marker_pos, mark_ypos1)
        ctx.line_to(marker_pos, mark_ypos2)
        ctx.set_source_rgba(0, 0, 0, 1)  # black
        ctx.set_line_width(8)
        ctx.stroke()
        # marker labels
        ctx.set_source_rgba(0, 0, 0, 1)
        ctx.set_font_size(55)
        ctx.move_to(marker_pos - 35, label_ypos)
        ctx.show_text(marker_label)
        ctx.stroke()


def make_rotated_labels(
    ctx, x, y, orientation, image_width, image_height, scaffold_name
):
    """Rotate the read names based on top/bottom orientation"""

    if orientation == "top":
        angle = -45
    elif orientation == "bottom":
        angle = 45

    # convert to radians
    theta = math.radians(angle)

    # get the mid point of the text if i were not rotated
    scaffold_name_midpoint = len(scaffold_name) / 2

    new_x = x + scaffold_name_midpoint + 25

    # if top or bottom chrom
    if angle > 0:
        new_y = y + scaffold_name_midpoint + 45
    else:
        new_y = y - scaffold_name_midpoint - 28

    # save un-rotated plot/context
    ctx.save()
    # create labels
    ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_source_rgb(0, 0, 0)  # black
    ctx.set_font_size(50)
    ctx.move_to(new_x, new_y)
    ctx.rotate(theta)
    ctx.show_text(scaffold_name)
    ctx.stroke()
    ctx.restore()


def Add_Sample_Names(cxt, sample_name, orientation, pdf_height, font_size=90):
    """Add the Names of the Focal Samples"""

    # establish top/bottom coordinates based
    if orientation == "top":
        x, y = 100, 75
    elif orientation == "bottom":
        x, y = 100, pdf_height - 50  # 50 positions above bottom of pdf

    cxt.set_source_rgb(0, 0, 0)  # black
    cxt.set_font_size(font_size)
    cxt.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    cxt.move_to(x, y)
    cxt.show_text(sample_name)
    cxt.stroke()


def Draw_Chrom_Borders(ctx, x_start, x_end, pos_y):
    """Draw borders around each scaffold"""

    # first horizontal line
    ctx.set_source_rgb(0, 0, 0)  # black
    ctx.move_to(x_start, pos_y + 23)
    ctx.line_to(x_end, pos_y + 23)
    ctx.set_line_width(5)
    ctx.stroke()

    # second horizontal line
    ctx.set_source_rgb(0, 0, 0)  # black
    ctx.move_to(x_start, pos_y - 23)
    ctx.line_to(x_end, pos_y - 23)
    ctx.set_line_width(5)
    ctx.stroke()

    # first veritcal line
    ctx.set_source_rgb(0, 0, 0)  # black
    ctx.move_to(x_start, pos_y + 23)
    ctx.line_to(x_start, pos_y - 23)
    ctx.set_line_width(5)
    ctx.stroke()

    # second vertical line
    ctx.set_source_rgb(0, 0, 0)  # black
    ctx.move_to(x_end, pos_y + 23)
    ctx.line_to(x_end, pos_y - 23)
    ctx.set_line_width(5)
    ctx.stroke()


def Draw_Chrome_Lines(
    ctx,
    sp_dict,
    scaffold_name,
    starting_position,
    genome_length,
    chrom_color,
    orientation,
    image_width,
):
    """Draw the chromosome lines"""

    # get length and estimate where it will end based on starting pos
    chrom_length = sp_dict[scaffold_name]
    chrom_end_point = (
        math.ceil((chrom_length / genome_length) * image_width) + starting_position
    )

    # draw the line
    ctx.move_to(starting_position, orientation)
    ctx.line_to(chrom_end_point, orientation)
    cairo_colors = quick_color_generator(chrom_color)
    ctx.set_source_rgba(
        cairo_colors.red,
        cairo_colors.green,
        cairo_colors.blue,
        cairo_colors.alpha,
    )
    ctx.set_line_width(50)
    ctx.stroke()

    # add border
    Draw_Chrom_Borders(ctx, starting_position, chrom_end_point, orientation)


def Draw_Ortho_Lines(
    ctx,
    ortho_row,
    starting_positionA,
    starting_positionB,
    genome_lengthA,
    genome_lengthB,
    coloring_iteration,
    top_chrom_pos,
    bottom_chrom_pos,
    image_width,
):
    """Draw The Line Between Orthologous Regions"""

    # parse ortho row list
    spA_scaf_length = ortho_row[0]
    spB_scaf_length = ortho_row[1]
    spA_start = ortho_row[2]
    spA_end = ortho_row[3]
    spB_start = ortho_row[4]
    spB_end = ortho_row[5]
    ortho_color = coloring_iteration % 8

    # get middle pos of each scaffold
    spA_pos = normalize_positions(
        starting_positionA,
        spA_start,
        spA_end,
        spA_scaf_length,
        genome_lengthA,
        image_width,
    )
    spB_pos = normalize_positions(
        starting_positionB,
        spB_start,
        spB_end,
        spB_scaf_length,
        genome_lengthB,
        image_width,
    )

    # P L O T
    ctx.move_to(spA_pos, top_chrom_pos + 25)
    ctx.line_to(spB_pos, bottom_chrom_pos - 25)
    cairo_colors = Pairing_Colors(ortho_color)
    ctx.set_source_rgba(
        cairo_colors.red,
        cairo_colors.green,
        cairo_colors.blue,
        cairo_colors.alpha,
    )
    ctx.set_line_width(2)
    ctx.stroke()


################ ATTENTION #################
def two_way_plot(
    spA,
    spB,
    full_dict,
    spA_dict,
    spB_dict,
    ortho_rows,
    color_scheme,
):
    """Create the two-way plot"""

    # create the image and context
    pdf_width = 7500
    pdf_height = 2000
    ims = cairo.PDFSurface("Testing.pdf", pdf_width, pdf_height)
    cairo_context = cairo.Context(ims)

    # constant variables (top & bottom chrom y axis values)
    y_top_chrom = 375
    y_bottom_chrom = 1625
    image_width = pdf_width - 400  # offset by starting positions on both sides

    # initiate x starting positions for both genomes
    current_starting_x = 100

    ######## Add sample names ########
    Add_Sample_Names(cairo_context, spA, "top", pdf_height)
    Add_Sample_Names(cairo_context, spB, "bottom", pdf_height)

    ###### calculate genome lengths #######
    spA_genome_length = sum(spA_dict.values())
    spB_genome_length = sum(spB_dict.values())

    ######## Get the starting and following chrom positions ########

    spA_starting_positions = establish_starting_positions(
        spA_dict, spA_genome_length, image_width, current_starting_x
    )

    spB_starting_positions = establish_starting_positions(
        spB_dict, spB_genome_length, image_width, current_starting_x
    )

    ######## Start looping through each chrom ########
    chroms_already_plottedA = []  # to avoid replotting chroms
    chroms_already_plottedB = []  # just in case scaffolds have same names

    ######## check color scheme #########
    if color_scheme == 1:
        current_iteration = False
    elif color_scheme == 2:
        current_iteration = 1
    else:
        sys.exit(f"Error: Color scheme must be either 1 or 2")

    # start plotting features
    for scaffold_names, ortho_rows in full_dict.items():
        spA_chrom = scaffold_names[0]
        spB_chrom = scaffold_names[1]
        current_starting_x_posA = spA_starting_positions[spA_chrom]
        current_starting_x_posB = spB_starting_positions[spB_chrom]
        # chrom colors
        spA_chrom_color = list(spA_starting_positions.keys()).index(spA_chrom)
        spB_chrom_color = list(spB_starting_positions.keys()).index(spB_chrom)

        # draw chrom lines
        if spA_chrom not in chroms_already_plottedA:

            Draw_Chrome_Lines(
                cairo_context,
                spA_dict,
                spA_chrom,
                current_starting_x_posA,
                spA_genome_length,
                spA_chrom_color,
                y_top_chrom,
                image_width,
            )

            chroms_already_plottedA.append(spA_chrom)

        if spB_chrom not in chroms_already_plottedB:
            Draw_Chrome_Lines(
                cairo_context,
                spB_dict,
                spB_chrom,
                current_starting_x_posB,
                spB_genome_length,
                spB_chrom_color,
                y_bottom_chrom,
                image_width,
            )
            chroms_already_plottedB.append(spB_chrom)

        ############## ADD scaffold names ###########

        make_rotated_labels(
            cairo_context,
            current_starting_x_posA,
            y_top_chrom,
            "top",
            image_width,
            pdf_height,
            spA_chrom,
        )

        make_rotated_labels(
            cairo_context,
            current_starting_x_posB,
            y_bottom_chrom,
            "bottom",
            image_width,
            pdf_height,
            spB_chrom,
        )

        # check for reverse
        reverse = check_for_reverse(ortho_rows)
        if reverse is True:
            ortho_rows = Reverse_Ortho_Rows(ortho_rows)

        # get color scheme
        if current_iteration is not False:
            coloring_iteration = current_iteration
            current_iteration += 1
        else:
            coloring_iteration = spA_chrom_color

        # add the pairing lines for that scaffold/chrom
        for ortho_row in ortho_rows:

            Draw_Ortho_Lines(
                cairo_context,
                ortho_row,
                current_starting_x_posA,
                current_starting_x_posB,
                spA_genome_length,
                spB_genome_length,
                coloring_iteration,
                y_top_chrom,
                y_bottom_chrom,
                image_width,
            )

    ######## Add tick marks ########
    Draw_Ticks(cairo_context, spA_genome_length, y_top_chrom, "top", image_width)
    Draw_Ticks(cairo_context, spB_genome_length, y_bottom_chrom, "bottom", image_width)

    # finish plotting
    cairo_context.show_page()


def three_way_plot(
    spA,
    spB,
    spC,
    full_dict_1v2,
    full_dict_2v3,
    spA_dict,
    spB_dict,
    spC_dict,
    ortho_rows_1v2,
    ortho_rows_2v3,
    color_scheme,
):
    # create the image and context
    pdf_width = 7500
    pdf_height = 4000
    ims = cairo.PDFSurface("Testing.pdf", pdf_width, pdf_height)
    cairo_context = cairo.Context(ims)

    # constant variables (top & bottom chrom y axis values)
    y_top_chrom = 375
    y_middle_chrom = 2000
    y_bottom_chrom = 3625
    image_width = pdf_width - 400  # offset by starting positions on both sides

    # all chroms will start 100 off
    current_starting_x = 100

    ######## Add sample names ########
    Add_Sample_Names(cairo_context, spA, "top", pdf_height)
    Add_Sample_Names(cairo_context, spC, "bottom", pdf_height)

    #### calculate genome lengths #####
    spA_genome_length = sum(spA_dict.values())
    spB_genome_length = sum(spB_dict.values())
    spC_genome_length = sum(spC_dict.values())

    ######## Get the starting and following chrom positions ########

    spA_starting_positions = establish_starting_positions(
        spA_dict, spA_genome_length, image_width, current_starting_x
    )

    spB_starting_positions = establish_starting_positions(
        spB_dict, spB_genome_length, image_width, current_starting_x
    )

    spC_starting_positions = establish_starting_positions(
        spC_dict, spC_genome_length, image_width, current_starting_x
    )

    ######## Start looping through each chrom ########
    chroms_already_plottedA = []  # to avoid replotting chroms
    chroms_already_plottedB = []  # just in case scaffolds have same names
    chroms_already_plottedC = []

    ######## check color scheme #########
    if color_scheme == 1:
        current_iteration = False
    elif color_scheme == 2:
        current_iteration = 1
    else:
        sys.exit(f"Error: Color scheme must be either 1 or 2")

    # start plotting spA vs spB
    for scaffold_names, ortho_rows in full_dict_1v2.items():
        spA_chrom = scaffold_names[0]
        spB_chrom = scaffold_names[1]
        current_starting_x_posA = spA_starting_positions[spA_chrom]
        current_starting_x_posB = spB_starting_positions[spB_chrom]
        # chrom colors
        spA_chrom_color = list(spA_starting_positions.keys()).index(spA_chrom)
        spB_chrom_color = list(spB_starting_positions.keys()).index(spB_chrom)

        # draw chrom lines
        if spA_chrom not in chroms_already_plottedA:

            Draw_Chrome_Lines(
                cairo_context,
                spA_dict,
                spA_chrom,
                current_starting_x_posA,
                spA_genome_length,
                spA_chrom_color,
                y_top_chrom,
                image_width,
            )

            chroms_already_plottedA.append(spA_chrom)

        if spB_chrom not in chroms_already_plottedB:
            Draw_Chrome_Lines(
                cairo_context,
                spB_dict,
                spB_chrom,
                current_starting_x_posB,
                spB_genome_length,
                spB_chrom_color,
                y_middle_chrom,
                image_width,
            )
            chroms_already_plottedB.append(spB_chrom)

        ############## ADD scaffold names for spA only ###########

        make_rotated_labels(
            cairo_context,
            current_starting_x_posA,
            y_top_chrom,
            "top",
            image_width,
            pdf_height,
            spA_chrom,
        )

        # check for reverse
        reverse = check_for_reverse(ortho_rows)
        if reverse is True:
            ortho_rows = Reverse_Ortho_Rows(ortho_rows)

        # get color scheme
        if current_iteration is not False:
            coloring_iteration = current_iteration
            current_iteration += 1
        else:
            coloring_iteration = spA_chrom_color

        # add the pairing lines for that scaffold/chrom
        for ortho_row in ortho_rows:

            Draw_Ortho_Lines(
                cairo_context,
                ortho_row,
                current_starting_x_posA,
                current_starting_x_posB,
                spA_genome_length,
                spB_genome_length,
                coloring_iteration,
                y_top_chrom,
                y_middle_chrom,
                image_width,
            )

    ######## Add tick marks for spA only ########
    Draw_Ticks(cairo_context, spA_genome_length, y_top_chrom, "top", image_width)

    # repeat with spB vs spC
    if color_scheme == 1:
        current_iteration = False
    elif color_scheme == 2:
        current_iteration = 1

    for scaffold_names, ortho_rows in full_dict_2v3.items():
        spB_chrom = scaffold_names[0]
        spC_chrom = scaffold_names[1]
        current_starting_x_posB = spB_starting_positions[spB_chrom]
        current_starting_x_posC = spC_starting_positions[spC_chrom]
        # chrom colors
        spB_chrom_color = list(spB_starting_positions.keys()).index(spB_chrom)
        spC_chrom_color = list(spC_starting_positions.keys()).index(spC_chrom)

        # draw chrom lines
        if spB_chrom not in chroms_already_plottedB:
            # will not plot spB chroms that are not present in spA vs spB comparison
            continue
        elif spC_chrom not in chroms_already_plottedC:
            Draw_Chrome_Lines(
                cairo_context,
                spC_dict,
                spC_chrom,
                current_starting_x_posC,
                spC_genome_length,
                spC_chrom_color,
                y_bottom_chrom,
                image_width,
            )
            chroms_already_plottedC.append(spC_chrom)

        ############## ADD scaffold names ###########

        make_rotated_labels(
            cairo_context,
            current_starting_x_posC,
            y_bottom_chrom,
            "bottom",
            image_width,
            pdf_height,
            spC_chrom,
        )

        # check for reverse
        reverse = check_for_reverse(ortho_rows)
        if reverse is True:
            ortho_rows = Reverse_Ortho_Rows(ortho_rows)

        # get color scheme
        if current_iteration is not False:
            coloring_iteration = current_iteration
            current_iteration += 1
        else:
            coloring_iteration = spB_chrom_color

        # add the pairing lines for that scaffold/chrom
        for ortho_row in ortho_rows:

            Draw_Ortho_Lines(
                cairo_context,
                ortho_row,
                current_starting_x_posB,
                current_starting_x_posC,
                spB_genome_length,
                spC_genome_length,
                coloring_iteration,
                y_middle_chrom,
                y_bottom_chrom,
                image_width,
            )

    ######## Add tick marks for bottom species ########
    Draw_Ticks(cairo_context, spC_genome_length, y_bottom_chrom, "bottom", image_width)

    # finish plotting
    cairo_context.show_page()


def run_two_way_plot(spA, spB, synolog_dir, color_scheme):
    """run the pipeline for a two way comparison"""

    # parse the membership file
    ortho_rows = get_ortho_rows(spA, spB, synolog_dir)
    # parse the chromosomes file
    spA_dict, spB_dict = get_scaffold_lengths(synolog_dir, spA, spB)
    # reorder chroms for spB
    spB_dict = chrom_reorder(spA_dict, spB_dict, ortho_rows)
    # merge objects to one big dictionary
    full_dict = merge_dicts(spA_dict, spB_dict, ortho_rows)
    # plot
    two_way_plot(
        spA,
        spB,
        full_dict,
        spA_dict,
        spB_dict,
        ortho_rows,
        color_scheme,
    )


def run_three_way_plot(spA, spB, spC, synolog_dir, color_scheme):
    """run the pipeline for a three way comparison"""

    # parse the membership files
    ortho_rows_1v2 = get_ortho_rows(spA, spB, synolog_dir)
    ortho_rows_2v3 = get_ortho_rows(spB, spC, synolog_dir)

    # parse the chromosome files
    spA_dict, spB_dict = get_scaffold_lengths(synolog_dir, spA, spB)
    # don't need redundant objects
    _, spC_dict = get_scaffold_lengths(synolog_dir, spB, spC)

    # re-orientate chrom dicts
    spB_dict = chrom_reorder(spA_dict, spB_dict, ortho_rows_1v2)
    spC_dict = chrom_reorder(spB_dict, spC_dict, ortho_rows_2v3)
    # merge objects to one big dictionary
    full_dict_1v2 = merge_dicts(spA_dict, spB_dict, ortho_rows_1v2)
    full_dict_2v3 = merge_dicts(spB_dict, spC_dict, ortho_rows_2v3)

    # plot three way comparison
    three_way_plot(
        spA,
        spB,
        spC,
        full_dict_1v2,
        full_dict_2v3,
        spA_dict,
        spB_dict,
        spC_dict,
        ortho_rows_1v2,
        ortho_rows_2v3,
        color_scheme,
    )


def main():
    """Put all the ingredients together and hope for the best!"""

    # get args
    spA, spB, spC, synolog_dir, color_scheme = get_arguments()

    # determine which pipeline to run
    if spC is None:
        run_two_way_plot(spA, spB, synolog_dir, color_scheme)
    else:
        run_three_way_plot(spA, spB, spC, synolog_dir, color_scheme)


if __name__ == "__main__":
    main()
