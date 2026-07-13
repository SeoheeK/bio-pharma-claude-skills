"""Command-line interface for the PK chatbot.

Usage:
  python -m pkchat.cli                    # interactive REPL
  python -m pkchat.cli "simulate CL=5 V1=30 dose=100 iv"   # one-shot
  python -m pkchat.cli --no-claude ...    # force the local router
"""
from __future__ import annotations

import argparse
import sys

from .chat.agent import PKChatAgent


def main(argv=None):
    ap = argparse.ArgumentParser(description="Chatbot-driven PK modelling engine")
    ap.add_argument("message", nargs="*", help="one-shot request (omit for REPL)")
    ap.add_argument("--no-claude", action="store_true",
                    help="force the local rule-based router")
    ap.add_argument("--model", default="claude-opus-4-8",
                    help="Claude model id (when the Claude backend is active)")
    args = ap.parse_args(argv)

    agent = PKChatAgent(model=args.model,
                        use_claude=False if args.no_claude else None)

    if args.message:
        print(agent.chat(" ".join(args.message)))
        return 0

    print(f"pkchat REPL  [backend: {agent.backend}]  (type 'quit' to exit)")
    print(agent.chat("help"))
    while True:
        try:
            line = input("\npk> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break
        if line.lower() in ("quit", "exit", "q"):
            break
        if not line:
            continue
        print(agent.chat(line))
    return 0


if __name__ == "__main__":
    sys.exit(main())
