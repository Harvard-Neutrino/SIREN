"""Load a CCM geometry product, a GDML that carries its own frame, into SIREN.

ccm_geometry writes every product with a <userinfo> block named ccm_product
that records, among other things, the detector's placement in target_sim
axes. A target_sim product is in those axes, so SIREN's geometry coordinates
are target_sim's and its detector coordinates are CCMAnalysis's once that
placement is applied. This module applies it.
"""
from .detector import DetectorModel, GeometryPosition
from .math import Matrix3D, Quaternion, Vector3D


def _find(auxiliaries, auxtype):
    found = [a for a in auxiliaries if a.type == auxtype]
    if len(found) != 1:
        raise ValueError(f"expected one {auxtype} auxiliary, found {len(found)}")
    return found[0]


def product_info(model):
    """The ccm_product block of the loaded GDML as a dict of its fields."""
    product = _find(model.GetGDMLUserinfo(), "ccm_product")
    fields = {a.type: a.value for a in product.children}
    placement = _find(product.children, "detector_to_target_sim")
    fields["detector_to_target_sim"] = {
        "status": placement.value,
        "rotation": [float(x) for x in _find(placement.children, "rotation").value.split()],
        "translation_m": [float(x) for x in _find(placement.children, "translation_m").value.split()],
    }
    fields["product"] = product.value
    return fields


def load_detector(path, strict=True):
    """A DetectorModel on the product at path, with the CCMAnalysis detector frame set."""
    model = DetectorModel()
    model.LoadGDML(str(path), strict)
    info = product_info(model)
    if info["frame"] != "target_sim":
        raise ValueError(f"SIREN loads products in target_sim axes; this one is in {info['frame']}")
    placement = info["detector_to_target_sim"]
    model.DetectorOrigin = GeometryPosition(Vector3D(*placement["translation_m"]))
    rotation = Quaternion()
    rotation.SetMatrix(Matrix3D(*placement["rotation"]))
    model.DetectorRotation = rotation
    return model
