import pdfplumber
import csv
import os

# Define function to search for tables with "adverse events" in caption
def extract_adverse_event_tables(pdf_path, output_dir):
    tables = []
    file_name = os.path.splitext(os.path.basename(pdf_path))[0]
    
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            text = page.extract_text()
            if text and ("adverse events" in text.lower() or "aes" in text.lower()):
                print(f"Found 'adverse events' on page {page.page_number} in {file_name}.")
                for table_num, table in enumerate(page.extract_tables(), start=1):
                    # Save each table to CSV
                    csv_filename = f"{output_dir}/{file_name}_table_{table_num}_page_{page.page_number}.csv"
                    with open(csv_filename, mode='w', newline='', encoding='utf-8') as csv_file:
                        writer = csv.writer(csv_file)
                        writer.writerows(table)
                    print(f"Exported table {table_num} from page {page.page_number} to {csv_filename}.")
                    tables.append(table)
                    
    return tables

# Paths to your PDF files
pdf_paths = [
    "38819031.pdf",
    "39282663.pdf"
]

# Directory to save CSV files
output_dir = "adverse_events_tables"
os.makedirs(output_dir, exist_ok=True)

# Extract tables with "adverse events" caption
for pdf_path in pdf_paths:
    tables = extract_adverse_event_tables(pdf_path, output_dir)
    if tables:
        print(f"Extracted and saved tables from {pdf_path}.")
    else:
        print(f"No 'adverse events' tables found in {pdf_path}.")
