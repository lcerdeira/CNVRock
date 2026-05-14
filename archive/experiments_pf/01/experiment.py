"""Thin wrapper — delegates to the central train.py with this experiment's config.

Preferred usage:  python ../../train.py config.yaml
                  (or via run_mac.sh / run_cluster.sh)
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "../.."))
sys.argv = [sys.argv[0], os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.yaml")]

from train import main
main()
