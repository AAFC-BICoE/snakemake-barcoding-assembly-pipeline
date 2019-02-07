import asyncio
import argparse
import csv
import json
from typing import Dict, Optional
from urllib.parse import urlencode

from Bio import SeqIO
import dataset
import requests

from engine import generate_output_content, parse_id_engine_xml, HEADERS


DATABASE_URL = "sqlite:///bold.sqlite"
DB = dataset.connect(DATABASE_URL)


def create_output_file(input_filename: str) -> str:
    """Containing only column headers of the CSV file."""
    output_filename = input_filename.strip() + "_output.csv"
    with open(output_filename, "w") as handle:
        csv_writer = csv.DictWriter(handle, fieldnames=HEADERS)
        csv_writer.writeheader()
    return output_filename


async def execute_one_record(db, output_filename, seq_record):
    print(f"* Reading seq {seq_record.name}")
    response = id_engine(seq_record, db, output_filename)
    seq_record_identifications = parse_id_engine_xml(response.text)
    # add our seq id to the list of identifications
    for seq_record_identification in seq_record_identifications:
        seq_record_identification["OtuID"] = seq_record.id
        taxonomy = get_taxonomy(seq_record_identification)
        try:
            bin = get_bin([seq_record_identification["ID"]])
        except KeyError:
            bin = ""
        seq_record_identification["BIN"] = bin
        seq_record_identification.update(taxonomy)
        seq_record_identification["seq_record"] = seq_record
    generate_output_content(seq_record_identifications, output_filename, seq_record)


def id_engine(seq_record, db, output_filename):
    """Send a COI sequence to BOLD and retrieve its identification"""
    print(f"* Processing sequence for {seq_record.id}")

    domain = "http://boldsystems.org/index.php/Ids_xml"
    payload = {'db': db, 'sequence': str(seq_record.seq)}
    url = domain + '?' + urlencode(payload)

    res = requests.get(url, headers={'User-Agent': 'bold_retriever'})
    return res


def get_taxonomy(seq_record: Dict[str, str]) -> Optional[Dict[str, str]]:
    tax_id = get_tax_id(seq_record)
    if tax_id:
        taxonomy = get_higher_level_taxonomy(tax_id)
        return taxonomy


def get_tax_id(seq_record: Dict[str, str]):
    tax_id = get_tax_id_from_db(seq_record)
    if tax_id:
        return tax_id

    domain = "http://boldsystems.org/index.php/API_Tax/TaxonSearch?taxName="
    url = domain + seq_record["taxonomicidentification"]
    res = requests.get(url, headers={'User-Agent': 'bold_retriever'})
    response_json = res.json()
    try:
        tax_id = response_json["top_matched_names"][0]["taxid"]
    except (KeyError, IndexError):
        tax_id = None

    if tax_id:
        table = DB["tax_ids"]
        data = {
            "taxon": seq_record["taxonomicidentification"],
            "tax_id": tax_id,
        }
        table.insert(data)
    return tax_id


def get_tax_id_from_db(seq_record: Dict[str, str]) -> Optional[str]:
    table = DB["tax_ids"]
    element = table.find_one(taxon=seq_record["taxonomicidentification"])
    if element:
        return element["tax_id"]


def get_higher_level_taxonomy(tax_id):
    table = DB["taxonomy"]
    element = table.find_one(tax_id=tax_id)
    if element:
        del element["id"]
        return element

    url = f"http://boldsystems.org/index.php/API_Tax/TaxonData?taxId={tax_id}" \
          "&dataTypes=basic&includeTree=true"
    res = requests.get(url, headers={'User-Agent': 'bold_retriever'})
    response_json = res.json()
    taxonomy = dict()

    for id in response_json.keys():
        category = response_json[id]
        value = category["taxon"]
        key = category["tax_rank"]
        taxonomy[key] = value

    taxonomy["tax_id"] = tax_id
    table.insert(taxonomy)
    return taxonomy


def get_bin(ids):
    ids = "|".join(ids)
    url = f"http://boldsystems.org/index.php/API_Public/specimen?ids={ids}&format=json"
    res = requests.get(url, headers={'User-Agent': 'bold_retriever'})

    records = None
    elements = None
    bin = None

    try:
        records = res.json()
    except json.decoder.JSONDecodeError:
        pass

    if records:
        try:
            elements = records["bold_records"]["records"]
        except KeyError:
            pass

    if elements:
        key = list(elements.keys())[0]

        try:
            bin = elements[key]["bin_uri"]
        except KeyError:
            pass

    return bin


async def main():
    if __name__ == "__main__":
        parser = argparse.ArgumentParser(description="bold_retriever")
        parser.add_argument('-f', '--filename', type=str, help='Fasta filename', required=True)
        parser.add_argument(
            '-db',
            '--database',
            choices=['COX1_SPECIES', 'COX1', 'COX1_SPECIES_PUBLIC', 'COX1_L640bp'],
            help='Choose a BOLD database. Enter one option.',
            required=True,
        )
        args = parser.parse_args()

        # Send seqs to BOLD Systems API and retrieve results
        fasta_file = args.filename
        db = args.database
        output_filename = create_output_file(fasta_file)
        print(f"Reading sequences from {fasta_file}")

        tasks = []
        for seq_record in SeqIO.parse(fasta_file, "fasta"):
            tasks.append(asyncio.ensure_future(execute_one_record(db, output_filename, seq_record)))
        await asyncio.gather(*tasks)


loop = asyncio.get_event_loop()
loop.run_until_complete(main())
loop.close()
