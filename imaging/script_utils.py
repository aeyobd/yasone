import tomllib


def read_img_keys():
    with open("image_keys.toml", "rb") as f:
        img_keys = tomllib.load(f)

    return img_keys



def calibration_folder(foldername):
    if foldername in ["20230708", "20230709", "20230723", "20230816"]:
        return foldername

    return {
        "yasone1": "20230708",
        "yasone2": "20230709",
        "yasone3a": "20230723",
        "yasone3b": "20230816",
    }[foldername]
