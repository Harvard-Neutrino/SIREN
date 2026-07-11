#!/usr/bin/env python3

"""Regression tests for exact directed-tiling volume normalization."""

import math
from types import SimpleNamespace
import unittest
from unittest import mock

import pytest
siren = pytest.importorskip("siren")
from siren import directed_tiling

class Vector3D:
    def __init__(self, x=0.0, y=0.0, z=0.0):
        self.x = x
        self.y = y
        self.z = z

    def GetX(self):
        return self.x

    def GetY(self):
        return self.y

    def GetZ(self):
        return self.z


class Quaternion:
    def __init__(self, x=0.0, y=0.0, z=0.0, w=1.0):
        self.X = x
        self.Y = y
        self.Z = z
        self.W = w


class Placement:
    def __init__(self, position=None, quaternion=None):
        self.Position = position or Vector3D()
        self.Quaternion = quaternion or Quaternion()


class Box:
    def __init__(self, *args):
        if len(args) == 3:
            self.placement = Placement()
            x, y, z = args
        elif len(args) == 4:
            self.placement, x, y, z = args
        else:
            raise TypeError("Box expects (x, y, z) or (placement, x, y, z)")
        self.X = x
        self.Y = y
        self.Z = z


class Cylinder:
    def __init__(self, radius, inner_radius, z, delta_phi):
        self.Radius = radius
        self.InnerRadius = inner_radius
        self.Z = z
        self.DeltaPhi = delta_phi


class Sphere:
    def __init__(self, radius, inner_radius, start_theta, delta_theta,
                 delta_phi):
        self.Radius = radius
        self.InnerRadius = inner_radius
        self.StartTheta = start_theta
        self.DeltaTheta = delta_theta
        self.DeltaPhi = delta_phi


class Composite:
    def __init__(self, volume=None):
        self.volume = volume


class Channel:
    def __init__(self, geometry, volume=None):
        self.geometry = geometry
        self.volume = volume


class MultiChannelPhaseSpace:
    def __init__(self):
        self.channels = []
        self.weights = []


class NestedMixtureChannel:
    def __init__(self, inner):
        self.inner = inner


class ConfigurationError(ValueError):
    pass


class DirectedTilingVolumeTest(unittest.TestCase):
    def setUp(self):
        def cpp_geometry_volume(geometry):
            if isinstance(geometry, Box):
                return geometry.X * geometry.Y * geometry.Z
            if isinstance(geometry, Cylinder):
                outer, inner = geometry.Radius, geometry.InnerRadius
                delta_phi = max(
                    0.0, min(2.0 * math.pi, geometry.DeltaPhi))
                return 0.5 * delta_phi * (
                    outer * outer - inner * inner) * geometry.Z
            if isinstance(geometry, Sphere):
                outer, inner = geometry.Radius, geometry.InnerRadius
                delta_phi = max(
                    0.0, min(2.0 * math.pi, geometry.DeltaPhi))
                theta_start = max(
                    0.0, min(math.pi, geometry.StartTheta))
                theta_end = max(theta_start, min(
                    math.pi, theta_start + geometry.DeltaTheta))
                return ((outer ** 3 - inner ** 3) / 3.0 * delta_phi
                        * (math.cos(theta_start) - math.cos(theta_end)))
            raise ValueError(
                "exact volume is unavailable for this geometry; provide volume_fn")

        self.cpp_geometry_volume = mock.Mock(side_effect=cpp_geometry_volume)
        fake_siren = SimpleNamespace(
            geometry=SimpleNamespace(
                Box=Box,
                Cylinder=Cylinder,
                Sphere=Sphere,
                Placement=Placement,
            ),
            math=SimpleNamespace(Vector3D=Vector3D),
            utilities=SimpleNamespace(ConfigurationError=ConfigurationError),
            injection=SimpleNamespace(
                geometry_volume=self.cpp_geometry_volume,
                MultiChannelPhaseSpace=MultiChannelPhaseSpace,
                NestedMixtureChannel=NestedMixtureChannel,
            ),
        )
        self.siren_patch = mock.patch.object(
            directed_tiling, "_siren", return_value=fake_siren)
        self.siren_patch.start()
        self.addCleanup(self.siren_patch.stop)

    def test_analytic_volumes_include_angular_cuts(self):
        box = Box(2.0, 3.0, 4.0)
        cylinder = Cylinder(5.0, 2.0, 8.0, math.pi / 3.0)
        sphere = Sphere(5.0, 2.0, math.pi / 4.0, math.pi / 3.0,
                        math.pi / 2.0)

        self.assertAlmostEqual(directed_tiling.geometry_volume(box), 24.0)
        self.assertAlmostEqual(
            directed_tiling.geometry_volume(cylinder),
            0.5 * (math.pi / 3.0) * (5.0 ** 2 - 2.0 ** 2) * 8.0,
        )
        self.assertAlmostEqual(
            directed_tiling.geometry_volume(sphere),
            ((5.0 ** 3 - 2.0 ** 3) / 3.0 * (math.pi / 2.0)
             * (math.cos(math.pi / 4.0)
                - math.cos(math.pi / 4.0 + math.pi / 3.0))),
        )
        self.assertEqual(self.cpp_geometry_volume.call_count, 3)

    def test_composite_volume_requires_exact_callback(self):
        with self.assertRaisesRegex(ValueError, "provide volume_fn"):
            directed_tiling.geometry_volume(Composite(), n_mc=1)

        with self.assertRaisesRegex(ValueError, "exact volume_fn"):
            directed_tiling._wrap(
                [Composite(2.0)], Channel, set_volume=True)

    def test_exact_callback_sets_volumes_and_drops_empty_tiles(self):
        tiles = [Composite(2.0), Composite(0.0), Composite(6.0)]

        wrapped = directed_tiling._wrap(
            tiles,
            Channel,
            set_volume=True,
            volume_fn=lambda geometry: geometry.volume,
        )

        self.assertIsInstance(wrapped, NestedMixtureChannel)
        self.assertEqual(
            [channel.geometry for channel in wrapped.inner.channels],
            [tiles[0], tiles[2]],
        )
        self.assertEqual(
            [channel.volume for channel in wrapped.inner.channels],
            [2.0, 6.0],
        )
        self.assertEqual(wrapped.inner.weights, [0.5, 0.5])

    def test_one_argument_factory_reports_configuration_error(self):
        def one_argument_factory(geometry):
            return Channel(geometry)

        with self.assertRaisesRegex(
                ConfigurationError,
                r"factory\(geometry, exact_volume\)"):
            directed_tiling._wrap(
                [Composite(2.0)],
                one_argument_factory,
                set_volume=True,
                volume_fn=lambda geometry: geometry.volume,
            )

    def test_two_argument_factory_type_errors_propagate(self):
        def failing_factory(geometry, exact_volume):
            raise TypeError("inner factory failure")

        with self.assertRaisesRegex(TypeError, "inner factory failure"):
            directed_tiling._wrap(
                [Composite(2.0)],
                failing_factory,
                set_volume=True,
                volume_fn=lambda geometry: geometry.volume,
            )

    def test_uninspectable_factory_is_called_directly(self):
        # Extension-type factories (pybind constructors) often have no
        # introspectable signature; the probe must assume the two-argument
        # form and call them rather than rejecting them as one-argument.
        class UninspectableFactory:
            @property
            def __signature__(self):
                raise ValueError("signature is not introspectable")

            def __call__(self, geometry, exact_volume=None):
                return Channel(geometry, exact_volume)

        wrapped = directed_tiling._wrap(
            [Composite(2.0)],
            UninspectableFactory(),
            set_volume=True,
            volume_fn=lambda geometry: geometry.volume,
        )

        self.assertEqual(wrapped.volume, 2.0)

    def test_all_empty_tiles_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "no non-empty tiles"):
            directed_tiling._wrap(
                [Composite(0.0), Composite(0.0)],
                Channel,
                set_volume=True,
                volume_fn=lambda geometry: geometry.volume,
            )

    def test_invalid_exact_callback_volume_is_rejected(self):
        for volume in (-1.0, math.nan, math.inf):
            with self.subTest(volume=volume):
                with self.assertRaisesRegex(
                        ValueError, "finite and non-negative"):
                    directed_tiling._wrap(
                        [Composite(volume)],
                        Channel,
                        set_volume=True,
                        volume_fn=lambda geometry: geometry.volume,
                    )

    def test_grid_cells_accepts_identity_rotation(self):
        box = Box(
            Placement(Vector3D(10.0, 20.0, 30.0), Quaternion(w=-1.0)),
            4.0, 2.0, 6.0)

        cells = directed_tiling.grid_cells(box, (2, 1, 1))

        self.assertEqual(len(cells), 2)
        self.assertEqual([(cell.X, cell.Y, cell.Z) for cell in cells],
                         [(2.0, 2.0, 6.0), (2.0, 2.0, 6.0)])
        self.assertEqual(
            [cell.placement.Position.GetX() for cell in cells],
            [9.0, 11.0])
        self.assertEqual(
            [cell.placement.Position.GetY() for cell in cells],
            [20.0, 20.0])
        self.assertEqual(
            [cell.placement.Position.GetZ() for cell in cells],
            [30.0, 30.0])

    def test_grid_cells_rejects_rotated_box(self):
        angle = math.pi / 2.0
        box = Box(
            Placement(
                Vector3D(),
                Quaternion(z=math.sin(angle / 2.0),
                           w=math.cos(angle / 2.0))),
            4.0, 2.0, 6.0)

        with self.assertRaisesRegex(
                ConfigurationError, "non-identity placement rotation"):
            directed_tiling.grid_cells(box, (2, 1, 1))


if __name__ == "__main__":
    unittest.main()
