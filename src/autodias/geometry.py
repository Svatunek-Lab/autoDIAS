from __future__ import annotations

import numpy as np


HARTREE_TO_KCAL = 627.509


def distance(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def angle(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> float:
    ba = (a - b) / np.linalg.norm(a - b)
    bc = (c - b) / np.linalg.norm(c - b)
    return float(np.degrees(np.arccos(np.clip(np.dot(ba, bc), -1.0, 1.0))))


def dihedral(p0: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    b0 = -1.0 * (p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2
    b1 /= np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return float(np.degrees(np.arctan2(y, x)))


def centroid(points: np.ndarray) -> np.ndarray:
    return points.mean(axis=0)


def rmsd_after_alignment(a: np.ndarray, b: np.ndarray) -> float:
    a_centered = a - centroid(a)
    b_centered = b - centroid(b)
    covariance = a_centered.T @ b_centered
    u, _, vh = np.linalg.svd(covariance)
    if (np.linalg.det(u) * np.linalg.det(vh)) < 0.0:
        u[:, -1] *= -1.0
    rotation = u @ vh
    aligned = a_centered @ rotation
    delta = aligned - b_centered
    return float(np.sqrt((delta * delta).sum() / len(a)))
