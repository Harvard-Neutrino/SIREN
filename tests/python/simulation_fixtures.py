"""Local detector geometry for facade tests; no resource lookup or download."""
from functools import lru_cache
from pathlib import Path
import siren


@lru_cache(maxsize=None)
def offline_detector(name='IceCube'):
    path = Path(__file__).resolve().parents[2] / 'resources/detectors' / name / (name + '-v1')
    model = siren.detector.DetectorModel()
    model.LoadMaterialModel(str(path / 'materials.dat'))
    model.LoadDetectorModel(str(path / 'densities.dat'))
    return model


def offline_fiducial():
    path = Path(__file__).resolve().parents[2] / 'resources/detectors/IceCube/IceCube-v1/densities.dat'
    lines = path.read_text().splitlines()
    fiducial = next(line for line in lines if line.startswith('fiducial '))
    detector = next(line for line in lines if line.startswith('detector '))
    return siren.detector.DetectorModel.ParseFiducialVolume(fiducial, detector)
