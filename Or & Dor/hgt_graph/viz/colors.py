import hashlib
import random
from ..constants import COLOR_BOUNDS


def stable_color_from_string(s: str) -> str:
    """
    Generate a stable color hex code from a given string.
    :param s: The input string.
    :return: A hex color code.
    """
    digest: bytes = hashlib.md5(s.encode("utf-8")).digest()
    seed: int = int.from_bytes(digest, byteorder="big")
    rnd: random.Random = random.Random(seed)
    r: int = rnd.randint(*COLOR_BOUNDS)
    g: int = rnd.randint(*COLOR_BOUNDS)
    b: int = rnd.randint(*COLOR_BOUNDS)
    return f"#{r:02x}{g:02x}{b:02x}"


def rescale_weight(w: float, w_min: float, w_max: float, min_w: float, max_w: float):
    """
    Rescale a single weight to a specified range [min_w, max_w].
    :param w: The weight to rescale.
    :param w_min: The minimum weight in the original range.
    :param w_max: The maximum weight in the original range.
    :param min_w: The minimum weight in the target range.
    :param max_w: The maximum weight in the target range.
    :return: The rescaled weight.
    """
    if w_min == w_max:
        return (min_w + max_w) / 2
    return min_w + (w - w_min) * (max_w - min_w) / (w_max - w_min)
