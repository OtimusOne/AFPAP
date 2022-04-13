'''
# File name: AFPAP_p2rank_gallery.py
# Description: Python script for generating P2Rank pocket gallery
# Author: Maghiar Octavian
# Date: 04-04-2022
'''
import argparse
import base64
import colorsys
import logging
import pathlib
import pandas as pd
from PIL import Image, ImageChops


def create_pocket_color_spectrum(nr_colors=1,  sat=0.6,  hue_min=0.6,  hue_max=1.2):
    '''
    Generate color spectrum for p2rank pocket visualization from HSV to RGB.

    Parameters:
        nr_colors: Number of colors.
        sat: Saturation.
        hue_min: Min hue.
        hue_max: Max hue.
    Returns:
        table of RGB colors
    '''
    color_table = [0]*nr_colors
    if nr_colors == 1:
        color_table[0] = create_pocket_color_spectrum(2, sat, hue_min, hue_max)[0]
    else:
        hue_range = hue_max - hue_min
        for i in range(nr_colors):
            hue = hue_min + hue_range * (i / nr_colors)
            if hue > 1:
                hue -= 1
            if i % 2 == 0:
                color_table[i] = colorsys.hsv_to_rgb(hue, sat,  0.9)
            else:
                color_table[i] = colorsys.hsv_to_rgb(hue, sat,  0.7)
    return color_table


def crop_image_whitespace(image, base_height=480):
    '''
    Crop image whitespace around borders.

    Parameters:
        image: Image object.
        base_height: Resized image height in px.
    Returns:
        Cropped image object.
    '''
    background_image = Image.new(image.mode, image.size, image.getpixel((0, 0)))
    diff = ImageChops.difference(image, background_image)
    diff = ImageChops.add(diff, diff, 2.0, -100)
    bbox = diff.getbbox()
    if bbox:
        image = image.crop(bbox)
    hpercent = (base_height/float(image.size[1]))
    wsize = int((float(image.size[0])*float(hpercent)))
    image = image.resize((wsize, base_height), Image.ANTIALIAS)
    return image


def main():
    '''
    Python script for generating P2Rank pocket viewer MultiQC report section.
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbosity', action="count", help="verbosity")
    parser.add_argument('-o', '--outputDir', type=pathlib.Path, default="./output", help="Output Directory")
    parser.add_argument('--AFPAPpath', type=pathlib.Path, required=True, help="Path to AFPAP home")

    args = parser.parse_args()
    console_logger = logging.StreamHandler()
    console_logger.setLevel(logging.INFO if args.verbosity else logging.ERROR)
    file_logger = logging.FileHandler(f"{args.outputDir}/workflow.log", mode='a')
    file_logger.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(name)s] %(levelname)-8.8s - %(message)s",
        handlers=[
            file_logger,
            console_logger
        ]
    )
    logging.info("P2Rank Pocket Viewer...")

    p2rank_predictions = pd.read_csv(f"{args.outputDir}/work/p2rank_predictions.csv")
    if p2rank_predictions.shape[0] == 0:
        return

    pocket_images = ["p2rank_1.png", "p2rank_2.png", "p2rank_3.png",
                     "p2rank_4.png", "p2rank_5.png", "p2rank_6.png",
                     "p2rank_7.png", "p2rank_8.png", "p2rank_9.png", "p2rank_10.png"]
    pocket_images_crop = ["p2rank_crop_1.png", "p2rank_crop_2.png", "p2rank_crop_3.png",
                          "p2rank_crop_4.png", "p2rank_crop_5.png", "p2rank_crop_6.png", "p2rank_crop_7.png", "p2rank_crop_8.png", "p2rank_crop_9.png", "p2rank_crop_10.png"]

    for img, img_crop in zip(pocket_images, pocket_images_crop):
        image_object = Image.open(f'{args.outputDir}/work/visualizations/{img}')
        image_object = crop_image_whitespace(image_object)
        image_object.save(f'{args.outputDir}/work/visualizations/{img_crop}', "PNG")

    with open(f"{args.AFPAPpath}/config/pocketViewer_template.html", 'r', encoding="utf8") as template:
        template_data = template.read()
        colors = create_pocket_color_spectrum(p2rank_predictions.shape[0], 0.6, 0.6, 1.2)
        pocket_description = "<span class=\"pocketDescription\">"
        for i in range(p2rank_predictions.shape[0]):
            pocket_description += f"<span class=\"pocketDescriptionSpan\">Pocket{i+1}<span class=\"pocketDescriptionBullet\" style=\"color:rgb({colors[i][0]*255},{colors[i][1]*255},{colors[i][2]*255});\">&#11044</span></span>\t "
        pocket_description += "</span>"
        template_data = template_data.replace("--pocketColors--", pocket_description)

        for j, img in enumerate(pocket_images_crop):
            with open(f'{args.outputDir}/work/visualizations/{img}', "rb") as image_file:
                img_data = base64.b64encode(image_file.read())
                template_data = template_data.replace(f"--pocket{j+1}--", img_data.decode('utf-8'))

        with open(f"{args.outputDir}/work/multiqc_files/p2rank_viewer_mqc.html", 'w', encoding="utf8") as mqc_file:
            mqc_file.write(template_data)


if __name__ == '__main__':
    main()
