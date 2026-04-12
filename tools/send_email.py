"""
Send a plain-text email via Gmail SMTP and print the Message-ID.

Usage:
    python tools/send_email.py --subject "..." --body "..." [--save-id path] [--in-reply-to "<msgid>"]

Reads EMAIL_ADDRESS and EMAIL_PASSWORD from .env at the repo root.
The Message-ID is printed to stdout so the caller can record it for
reply tracking (check_reply.py).

Pass --in-reply-to to thread the email under an existing conversation.
"""

import argparse
import os
import smtplib
import sys
from email.message import EmailMessage
from email.utils import make_msgid


def load_env(repo_root):
    env = {}
    env_path = os.path.join(repo_root, ".env")
    with open(env_path) as f:
        for line in f:
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env


def send(subject, body, save_id_path=None, in_reply_to=None):
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env = load_env(repo_root)
    address  = env["EMAIL_ADDRESS"]
    password = env["EMAIL_PASSWORD"]

    msg = EmailMessage()
    msg["From"]           = address
    msg["To"]             = address
    msg["Subject"]        = subject
    msg["Message-ID"]     = make_msgid(domain="deep-gw-cnv")
    msg["X-CNV-Daemon"]   = "ack"   # marks all our outgoing emails so the reply checker ignores them
    if in_reply_to:
        msg["In-Reply-To"] = in_reply_to
        msg["References"]  = in_reply_to  # ensures Gmail threads correctly
    msg.set_content(body)

    with smtplib.SMTP("smtp.gmail.com", 587) as smtp:
        smtp.starttls()
        smtp.login(address, password)
        smtp.send_message(msg)

    msg_id = msg["Message-ID"]
    print(f"Sent. Message-ID: {msg_id}")

    if save_id_path:
        with open(save_id_path, "w") as f:
            f.write(msg_id)
        print(f"Message-ID saved to {save_id_path}")

    return msg_id


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--subject",     required=True)
    parser.add_argument("--body",        required=True, help="Body text or @path to read from file")
    parser.add_argument("--save-id",     default=None,  help="File to save the Message-ID into")
    parser.add_argument("--in-reply-to", default=None,  help="Message-ID to thread this email under")
    args = parser.parse_args()

    body = args.body
    if body.startswith("@"):
        with open(body[1:]) as f:
            body = f.read()

    send(args.subject, body, args.save_id, args.in_reply_to)
