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


def _decode_body(raw_bytes):
    """Parse an RFC822 message and return its decoded plain-text body, or ''."""
    msg = email.message_from_bytes(raw_bytes)
    for part in msg.walk():
        if part.get_content_type() == "text/plain":
            payload = part.get_payload(decode=True)  # handles base64/quoted-printable
            if payload:
                return payload.decode(part.get_content_charset() or "utf-8", errors="replace").strip()
    return ""


# Body prefixes that identify daemon-sent ack emails (for legacy emails that
# predate the X-CNV-Daemon header). These are matched case-insensitively
# against the start of the decoded body.
_ACK_BODY_PREFIXES = (
    "got it. i'll review your feedback",
    "experiment ",           # "Experiment N has started…" / "Experiment N has finished…"
    "automatic proposal failed",
)


def _is_own_ack(parsed, body):
    """Return True if this email is a daemon-sent ack we should skip."""
    if parsed.get("X-CNV-Daemon"):
        return True  # header present — definitive
    # Fallback for old acks sent before the header was added.
    if body:
        lower = body.lower()
        if any(lower.startswith(p) for p in _ACK_BODY_PREFIXES):
            return True
    return False


def _fetch_reply_body(imap, original_msg_id, own_address):
    """Search for and return the decoded plain-text body of the user's reply, or None.

    Iterates oldest-first. Skips any email that looks like a daemon ack
    (either via the X-CNV-Daemon header or known ack body prefixes).
    """
    _, data = imap.search(None, f'HEADER "In-Reply-To" "{original_msg_id}"')
    uids = data[0].split()
    if not uids:
        # Some clients set References instead of In-Reply-To
        _, data = imap.search(None, f'HEADER "References" "{original_msg_id}"')
        uids = data[0].split()
    if not uids:
        return None

    for uid in uids:
        _, msg_data = imap.fetch(uid, "(RFC822)")
        raw = msg_data[0][1]
        parsed = email.message_from_bytes(raw)
        body = _decode_body(raw)
        if _is_own_ack(parsed, body):
            continue  # daemon ack — skip
        return body

    return None  # only daemon acks found — treat as no reply


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
        body = _fetch_reply_body(imap, original_msg_id, address)

    if body is None:
        print("NO_REPLY")
        sys.exit(1)

    if print_body:
        # Used by /propose-experiment to read feedback. Body goes to stdout only.
        print(body)
        return

    # Body is intentionally NOT logged or written to disk.
    # Strip quoted lines ("> ...") so a forwarded/replied AUTHORISE in the
    # original proposal doesn't trigger a false positive.
    unquoted = "\n".join(
        line for line in body.splitlines() if not line.startswith(">")
    )
    if "AUTHORISE" in unquoted.upper():
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
