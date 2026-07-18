#!/usr/bin/env python3
"""Create a compact contact sheet from E2E preview images."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from PIL import Image, ImageDraw, ImageOps


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--title", default="")
    args = parser.parse_args()

    paths = sorted(path for path in args.input_dir.rglob("*.png") if path.resolve() != args.output.resolve())
    if not paths:
        parser.error(f"no PNG images found under {args.input_dir}")

    tile_width, image_height, caption_height = 760, 520, 44
    rows = math.ceil(len(paths) / args.columns)
    title_height = 56 if args.title else 0
    sheet = Image.new(
        "RGB",
        (tile_width * args.columns, title_height + (image_height + caption_height) * rows),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    if args.title:
        draw.text((16, 18), args.title, fill="black")

    for index, path in enumerate(paths):
        with Image.open(path) as source:
            image = ImageOps.contain(source.convert("RGB"), (tile_width - 24, image_height - 24))
        column, row = index % args.columns, index // args.columns
        x = column * tile_width + (tile_width - image.width) // 2
        y = title_height + row * (image_height + caption_height) + (image_height - image.height) // 2
        sheet.paste(image, (x, y))
        caption = path.relative_to(args.input_dir).as_posix()
        draw.text(
            (column * tile_width + 12,
             title_height + row * (image_height + caption_height) + image_height + 10),
            caption,
            fill="black",
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(args.output, dpi=(150, 150), optimize=True)
    print(f"Created {args.output} with {len(paths)} images")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
