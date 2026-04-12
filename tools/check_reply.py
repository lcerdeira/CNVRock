"""
Check Gmail inbox for a reply to one specific proposal email.

This script is deliberately narrow:
  - It searches ONLY by the exact Message-ID of the proposal email.
  - It fetches ONLY the single matching reply, and ONLY its body.
  - It never lists, reads, or logs any other email.
  - The reply body is never written to disk or to the daemon log.

Usage:
    python tools/check_reply.py --msg-id-file path/to/msgid.txt

Exits with:
    0  and prints "AUTHORISE"          — reply found containing "AUTHORISE"
    0  and prints "FEEDBACK"           — reply found with other content
    1  and prints "NO_REPLY"           — no matching reply yet

Reads EMAIL_ADDRESS / EMAIL_PASSWORD from .env at the repo root.
"""

import argparse
import email
import imaplib
import os
import sys


def load_env(repo_root):
    env = {}
    with open(os.path.join(repo_root, ".env")) as f:
        for line in f:
            line = line.strip()
            if "=" in line and not line.startswith("#"):
                k, v = line.split("=", 1)
                env[k.strip()] = v.strip()
    return env


def _fetch_reply_body(imap, original_msg_id):
    """Search for and return the decoded plain-text body of the reply, or None."""
    _, data = imap.search(None, f'HEADER "In-Reply-To" "{original_msg_id}"')
    uids = data[0].split()
    if not uids:
        # Some clients use References instead of In-Reply-To
        _, data = imap.search(None, f'HEADER "References" "{original_msg_id}"')
        uids = data[0].split()
    if not uids:
        return None
    # Fetch the full RFC822 message so we can decode MIME properly (e.g. base64).
    _, msg_data = imap.fetch(uids[-1], "(RFC822)")
    raw = msg_data[0][1]
    msg = email.message_from_bytes(raw)
    # Walk parts and return the first text/plain payload.
    for part in msg.walk():
        if part.get_content_type() == "text/plain":
            payload = part.get_payload(decode=True)  # handles base64/quoted-printable
            if payload:
                return payload.decode(part.get_content_charset() or "utf-8", errors="replace").strip()
    return ""


def check(msg_id_file, print_body=False):
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    env      = load_env(repo_root)
    address  = env["EMAIL_ADDRESS"]
    password = env["EMAIL_PASSWORD"]

    with open(msg_id_file) as f:
        original_msg_id = f.read().strip()

    with imaplib.IMAP4_SSL("imap.gmail.com") as imap:
        imap.login(address, password)
        imap.select("INBOX")
        body = _fetch_reply_body(imap, original_msg_id)

    if body is None:
        print("NO_REPLY")
        sys.exit(1)

    if print_body:
        # Used by /propose-experiment to read feedback. Body goes to stdout only.
        print(body)
        return

    # Body is intentionally NOT logged or written to disk.
    if "AUTHORISE" in body.upper():
        print("AUTHORISE")
    else:
        print("FEEDBACK")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--msg-id-file", required=True,
                        help="File containing the Message-ID of the proposal email")
    parser.add_argument("--print-body", action="store_true",
                        help="Print the reply body to stdout (for /propose-experiment feedback loop)")
    args = parser.parse_args()
    check(args.msg_id_file, print_body=args.print_body)
