from PyPDF2 import PdfFileWriter, PdfFileReader
import os

inputpdf = PdfFileReader(open(snakemake.input[0], "rb"))
out_dir = snakemake.output[0]
os.makedirs(out_dir, exist_ok=True)

for i in range(inputpdf.numPages):
    output = PdfFileWriter()
    output.addPage(inputpdf.getPage(i))
    with open(f"{snakemake.output[0]}/qc_{i:02d}.pdf", "wb") as outputStream:
        output.write(outputStream)
