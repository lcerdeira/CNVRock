#!/usr/bin/env python3
import os
import time
import smtplib
import subprocess
import argparse

from dotenv import load_dotenv
from email.message import EmailMessage
from datetime import datetime

def get_bpeek_tail(job_id, n=50) -> str:
    """Needs to be done as malpipe"""
    result = subprocess.run(
        ["bpeek", str(job_id)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if result.returncode != 0:
        return f"bpeek failed:\n{result.stderr}"

    lines = result.stdout.splitlines()
    return "\n".join(lines[-n:])

def check_job_is_alive(job_id) -> bool:
    cmd = ["bjobs", str(job_id)]
    result = subprocess.run(
        cmd, text=True, check=False,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # Log the result for debugging
    print(f"bjobs return code: {result.returncode}")
    print(f"bjobs stdout: {result.stdout}")
    print(f"bjobs stderr: {result.stderr}")

    if result.returncode != 0:
        print("Return code is not 0")
        return False

    if "EXIT" in result.stdout or "EXIT" in result.stderr:
        print("EXIT detected")
        return False

    if "Command error" in get_bpeek_tail(job_id):
        print("command error found in bpeek")
        return False

    return str(job_id) in result.stdout

def send_email(sender, password, subject, body):
    msg = EmailMessage()
    msg["From"] = sender
    msg["To"] = sender
    msg["Subject"] = subject
    msg.set_content(body)

    # Create fresh connection each time
    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(sender, password)
        server.send_message(msg)

def monitor_job(job_id, sender, password):
    job_was_alive = True
    last_daily_sent = None

    # Remove the unused SMTP connection here
    now = datetime.now()
    print(f"Monitoring started {now}.")

    send_email(
        sender, password,
        subject=f"LSF job {job_id} tracking begun {now}.",
        body=get_bpeek_tail(job_id),
    )

    while True:
        # Update 'now' in each loop iteration
        now = datetime.now()
        
        # ---- hourly job check ----
        is_alive = check_job_is_alive(job_id)

        if job_was_alive and not is_alive:
            send_email(
                sender, password,
                subject=f"LSF job {job_id} finished",
                body=get_bpeek_tail(job_id),
            )
            print("Job finished, alert sent.")
            return  # stop monitoring

        job_was_alive = is_alive

        # ---- daily 12:00 bpeek ----
        if now.hour == 6:
            if last_daily_sent != now.date():
                send_email(
                    sender, password,
                    subject=f"Daily bpeek update for job {job_id}",
                    body=get_bpeek_tail(job_id),
                )

                last_daily_sent = now.date()
                print("Daily update sent.")

        # sleep until next minute
        time.sleep(60)

def main():
    load_dotenv()

    parser = argparse.ArgumentParser(description="Monitor LSF job and send email alerts")
    parser.add_argument("job_id", help="LSF JOBID to track (USER must be malpipe)")
    args = parser.parse_args()

    """Create a Gmail App password at: https://myaccount.google.com/apppasswords"""
    sender = os.environ.get("EMAIL_ADDRESS")
    password = os.environ.get("EMAIL_PASSWORD")

    if not sender or not password:
        raise RuntimeError("Missing email credentials. Set EMAIL_ADDRESS and EMAIL_PASSWORD")

    job_id = args.job_id

    if not check_job_is_alive(job_id):
        raise RuntimeError(f"Job {job_id} not found or not running")

    monitor_job(job_id, sender, password)

if __name__ == "__main__":
    main()
