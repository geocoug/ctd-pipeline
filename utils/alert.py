import argparse
import os
import smtplib
from email.mime.text import MIMEText

import cf_xarray  # noqa
import xarray as xr
from dotenv import load_dotenv

load_dotenv()
PASSWORD = os.getenv("PASSWORD")
SUBJECT = "ASV CTD QC Alert"
SENDER = "grantcaleb22@gmail.com"
RECIPIENTS = ["grantcaleb22@gmail.com"]


def send_alert(
    subject: str,
    body: str,
    sender: str,
    recipients: list | tuple,
    password: str,
) -> None:
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = sender
    msg["To"] = ", ".join(recipients)
    try:
        smtp_server = smtplib.SMTP_SSL("smtp.gmail.com", 465)
        smtp_server.login(sender, password)
        smtp_server.sendmail(sender, recipients, msg.as_string())
        smtp_server.quit()
    except Exception:
        raise


def check_ncfile(input_file: str, ncfile: str) -> tuple[bool, str | None]:
    """Review output NetCDF QC file and alert critical staff of issues.

    Check that an ASCII file and NetCDF file exist
    Check that the ASCII and NetCDF have the same number of records?
        May need to change if file headers change.
    Check that XX parameter has failed?

    Args:
        input_file (str): _description_
        ncfile (str): _description_

    Returns:
        tuple[bool, str | None]: _description_
    """
    if not os.path.exists(ncfile):
        return True, f"NetCDF file not created for {input_file}."
    try:
        with open(input_file) as f:
            source_records = len(f.readlines()) - 1
    except Exception:
        return True, f"Unable to read observation file {input_file}"
    try:
        qc = xr.open_dataset(ncfile, decode_times=False)
        qc_records = qc.dims["time"]
    except Exception:
        return (
            True,
            f"NetCDF {ncfile} file is corrupt. Observations are in {input_file}.",
        )
    if source_records != qc_records:
        return (
            True,
            f"The number of records in the source file {input_file} (n={source_records}) and QC NetCDF {ncfile} (n={qc_records}) do not match.",  # noqa
        )
    return False, None


def main(input_file: str) -> None:
    if not os.path.exists(input_file):
        send_alert(
            SUBJECT,
            f"Observation file {input_file} is corrupt.",
            SENDER,
            RECIPIENTS,
            PASSWORD,
        )
    ncfile = f"{input_file}.nc".replace("received", "processed")
    error, message = check_ncfile(input_file, ncfile)
    if error:
        if message:
            send_alert(SUBJECT, message, SENDER, RECIPIENTS, PASSWORD)
    else:
        subject = "ASV CTD Cast Received"
        message = f"Cast observations for {input_file} have been received and a QC file has been generated at {ncfile}."  # noqa
        send_alert(subject, message, SENDER, RECIPIENTS, PASSWORD)


def clparser() -> argparse.ArgumentParser:
    """Create a parser to handle input arguments and displaying.

    a script specific help message.
    """
    desc_msg = """Alert critical staff of successful data exchange or issues."""
    parser = argparse.ArgumentParser(description=desc_msg)
    parser.add_argument(
        "input_file",
        help="Path to the input sensor data file.",
    )
    return parser


if __name__ == "__main__":
    parser = clparser()
    args = parser.parse_args()
    input_file = args.input_file
    main(input_file=input_file)
