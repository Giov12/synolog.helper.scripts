#!/bin/env python3

import argparse
import os
import pandas as pd

excel_f = '' # excel file
odir    = '.'

def get_arguments() -> int:
    """get the arguments"""

    global excel_f, odir

    d = "Write a csv file for each sheet in a comma delimited excel file"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-e", "--excel", help="excel file [required]", required=True, type=str)
    parser.add_argument("-o", "--outdir", help="directory to place the csv files", type=str, default='./')

    args    = parser.parse_args()
    excel_f = args.excel
    odir    = args.outdir

    assert os.path.isfile(excel_f), f"Could not locate file: {excel_f}"
    assert os.path.isdir(odir), f"Could not locate directory: {odir}"

    if (odir[-1] == '/'):
        odir = odir[:-1]

    return 0


def parse_excel() -> int:
    """this function will handle all the writing"""

    global odir, excel_f

    sheets = pd.read_excel(excel_f, sheet_name=None) # create a mapping of sheet -> df

    for sheet_name, df in sheets.items():

        # drop empty rows
        df  = df.dropna(how="all")
        out = f"{odir}/{sheet_name}.csv"
        df.to_csv(out, index = False)
        print("wrote", out)

    return 0

def main() -> int:
    """entry point to this tiny application"""

    # get arguments
    get_arguments()

    # parse & write
    parse_excel()

    return 0

if __name__ == "__main__":
    main()
