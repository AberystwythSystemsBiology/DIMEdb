import os
import logging

logging.basicConfig(level=logging.DEBUG)

from refs import *

data_directory = os.path.join(os.path.expanduser("~"), "Data/dimedb/")

hmdb_path = os.path.join(data_directory, "hmdb/")
hmdb_fp = os.path.join(hmdb_path, "hmdb_metabolites.zip")

massbank_path = os.path.join(data_directory, "massbank/")

spektraris_path = os.path.join(data_directory, "spektraris/")
spektraris_fp = os.path.join(spektraris_path, "database_flatfiles.zip")

nmrshiftdb_path = os.path.join(data_directory, "nmrshiftdb/")
nmrshiftdb_fp = os.path.join(nmrshiftdb_path, "nmrshiftdb2.sd")

respect_path = os.path.join(data_directory, "respect/")
respect_fp = os.path.join(respect_path, "respect.zip")
respect_udir = os.path.join(respect_path, "unzipped/")

golm_path = os.path.join(data_directory, "golm/")


if __name__ == "__main__":
    if os.path.isdir(data_directory) != True:
        os.makedirs(data_directory)

    # GOLM
    if os.path.isdir(golm_path) != True:
        os.makedirs(golm_path)
        golm.download(golm_path)
    if os.path.isfile(os.path.join(data_directory, "golm.json")) != True:
        golm.convert(golm_path, "golm.json")
    # Respect
    if os.path.isdir(respect_path) != True:
        os.makedirs(respect_path)
    if os.path.isfile(respect_fp) != True:
        respect.download(respect_fp)
    if os.path.isdir(respect_udir) != True:
        os.makedirs(respect_udir)
        respect.unzip(respect_fp, respect_udir)
    if os.path.isfile(os.path.join(data_directory, "respect.json")) != True:
        respect.process(respect_udir,  os.path.join(data_directory, "respect.json"))

    # NMR ShiftDB 2
    if os.path.isdir(nmrshiftdb_path) != True:
        os.makedirs(nmrshiftdb_path)
    if os.path.isfile(nmrshiftdb_fp) != True:
        nmrshiftdb.download(nmrshiftdb_fp, nmrshiftdb_path)
    if os.path.isfile(os.path.join(data_directory, "nmrshiftdb2.json")) != None:
        nmrshiftdb.convert(nmrshiftdb_fp, os.path.join(data_directory, "nmrshiftdb2.json"))

    # Spektraris
    if os.path.isdir(spektraris_path) != True:
        os.makedirs(spektraris_path)
    if os.path.isfile(spektraris_path) != True:
        spektraris.download(spektraris_path, spektraris_fp)
    spektraris.unzip(spektraris_path, spektraris_fp)
    if os.path.isfile(os.path.join(data_directory, "spektraris.json")) != True:
        spektraris.convert(os.path.join(spektraris_path, "record.tsv"), os.path.join(data_directory, "spektraris.json"))
    # Massbank
    if os.path.isdir(massbank_path) != True:
        os.makedirs(massbank_path)
    if os.path.isdir(os.path.join(massbank_path, "record")) != True:
        massbank.download(massbank_path)
    if os.path.isfile(os.path.join(data_directory, "massbank.json")) != True:
        massbank.convert_data(os.path.join(massbank_path, "record"), os.path.join(data_directory, "massbank.json"))
    # HMDB
    hmdb.download(hmdb_path, hmdb_fp)
    if os.path.isfile(os.path.join(hmdb_path, "hmdb_metabolites.xml")) != True:
        hmdb.unzip(hmdb_path, hmdb_fp)
    if os.path.isfile(os.path.join(data_directory, "hmdb.json")) != True:
        hmdb.extract_from_xml(os.path.join(hmdb_path, "hmdb_metabolites.xml"), os.path.join(data_directory, "hmdb.json"))
