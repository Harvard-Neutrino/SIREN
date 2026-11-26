"""Side panel of the interactive viewer: sections, layers and materials."""
from ._visualization_presets import ACCENT, INK, OTHER, PANEL, RULE, SOFT, clean_name

WIDTH = 270     # panel width in points; scaled by the display's pixel density
_PAD, _GAP = 14, 5
# Colors offered for a row; "Default" returns to the settings or density color.
PALETTE = [(.753, .278, .227), (.878, .569, .227), (.890, .765, .294), (.561, .769, .420),
           (.184, .612, .471), (.247, .608, .820), (.294, .384, .788), (.482, .408, .722),
           (.659, .345, .604), (.851, .478, .608), (.541, .353, .235), (.796, .741, .549),
           (.655, .678, .698), (.349, .376, .420), (.910, .875, .816), (.110, .165, .200)]


class Panel:
    """Draws rows in its own viewport and maps clicks back to actions.

    The session supplies the state; the panel only displays it. Rows are laid
    out from the top; the mouse wheel scrolls them when they do not fit.
    """

    def __init__(self, vtk, window, session):
        self.vtk, self.window, self.session = vtk, window, session
        self.ren = vtk.vtkRenderer()
        self.ren.SetBackground(*PANEL)
        self.ren.InteractiveOff()
        window.AddRenderer(self.ren)
        self.scroll, self.height = 0., 0.
        self.visible = True
        self.targets = []  # (x0, y0, x1, y1, action), window pixels
        self.drag = None
        self._props = []
        self._fonts = {}

    # ---- geometry -------------------------------------------------------
    @property
    def scale(self):
        return max(self.window.GetDPI() / 72., 1.)

    def pixel_width(self):
        if not self.visible:
            return 0
        return int(min(WIDTH * self.scale, .45 * self.window.GetSize()[0]))

    def place(self, *others):
        """Give the panel its strip of the window and the rest to *others*."""
        width = max(self.window.GetSize()[0], 1)
        edge = self.pixel_width() / width
        self.ren.SetDraw(self.visible)
        self.ren.SetViewport(0., 0., max(edge, 1e-6), 1.)
        for ren in others:
            ren.SetViewport(edge, 0., 1., 1.)

    def contains(self, x):
        return self.visible and x < self.pixel_width()

    # ---- drawing --------------------------------------------------------
    def _rect(self, x0, y0, x1, y1, color, filled=True):
        vtk = self.vtk
        points, cells = vtk.vtkPoints(), vtk.vtkCellArray()
        for x, y in ((x0, y0), (x1, y0), (x1, y1), (x0, y1)):
            points.InsertNextPoint(x, y, 0)
        order = (0, 1, 2, 3) if filled else (0, 1, 2, 3, 0)
        cells.InsertNextCell(len(order))
        for i in order:
            cells.InsertCellPoint(i)
        data = vtk.vtkPolyData()
        data.SetPoints(points)
        data.SetPolys(cells) if filled else data.SetLines(cells)
        mapper = vtk.vtkPolyDataMapper2D()
        mapper.SetInputData(data)
        actor = vtk.vtkActor2D()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(*color)
        self._add(actor)

    def _add(self, prop):
        self.ren.AddViewProp(prop)
        self._props.append(prop)

    def _text(self, text, x, y, size=12, color=INK, bold=False, center=False, right=False):
        key = (size, color, bold, center, right)
        if key not in self._fonts:
            style = self.vtk.vtkTextProperty()
            style.SetFontFamilyToArial()
            style.SetFontSize(size)
            style.SetColor(*color)
            style.SetBold(bold)
            style.SetVerticalJustificationToCentered()
            if center:
                style.SetJustificationToCentered()
            elif right:
                style.SetJustificationToRight()
            self._fonts[key] = style
        actor = self.vtk.vtkTextActor()
        actor.SetInput(text)
        actor.SetTextProperty(self._fonts[key])
        actor.SetPosition(x, y)
        self._add(actor)

    def _fit(self, text, width, size):
        # Arial averages a little over half its size per character.
        room = max(int(width / (.52 * size * self.scale)), 4)
        return text if len(text) <= room else text[:room - 3].rstrip() + "..."

    def _wrap(self, text, width, size):
        room, lines = max(int(width / (.5 * size * self.scale)), 8), [""]
        for word in text.split():
            if lines[-1] and len(lines[-1]) + 1 + len(word) > room:
                lines.append("")
            lines[-1] = (lines[-1] + " " + word).strip()
        return lines

    def _button(self, label, x0, y0, x1, y1, pressed, action):
        self._rect(x0, y0, x1, y1, ACCENT if pressed else PANEL)
        if not pressed:
            self._rect(x0, y0, x1, y1, RULE, filled=False)
        self._text(self._fit(label, x1 - x0 - 6, 11), (x0 + x1) / 2, (y0 + y1) / 2, 11,
                   (1., 1., 1.) if pressed else INK, center=True)
        self.targets.append((x0, y0, x1, y1, action))

    def rebuild(self):
        session, s = self.session, self.scale
        for prop in self._props:
            self.ren.RemoveViewProp(prop)
        self._props, self.targets = [], []
        if not self.visible:
            return
        width, window_height = self.pixel_width(), self.window.GetSize()[1]
        left, right = _PAD * s, width - _PAD * s
        self._rect(width - 1, 0, width, window_height, RULE)
        y = [window_height - 14 * s + self.scroll]

        def down(points):
            y[0] -= points * s
            return y[0]

        def heading(text):
            down(12)
            self._text(text, left, down(9), 14, bold=True)
            down(12)

        def buttons(rows, columns):
            gap, height = _GAP * s, 23 * s
            cell = (right - left - gap * (columns - 1)) / columns
            x, used = left, 0
            for label, pressed, action, span in rows:
                span = min(span, columns)
                if used + span > columns:
                    x, used = left, 0
                    y[0] -= height + gap
                w = cell * span + gap * (span - 1)
                self._button(label, x, y[0] - height, x + w, y[0], pressed, action)
                x, used = x + w + gap, used + span
            y[0] -= height

        def check(label, on, action, swatch=None, clear=False, kind=None, toggle=None):
            """A row: shown box, color swatch (hollow if none), name, see-through mark.

            The box always shows or hides (*toggle*); the rest of the row runs *action*.
            """
            middle, box = down(10), 6 * s
            self._rect(left, middle - box, left + 2 * box, middle + box, ACCENT if on else PANEL)
            if on:
                self._tick(left, middle, box)
            else:
                self._rect(left, middle - box, left + 2 * box, middle + box, SOFT, filled=False)
            start = left + 20 * s
            if kind is not None:
                if swatch is not None:
                    self._rect(start, middle - box, start + 2 * box, middle + box, swatch)
                self._rect(start, middle - box, start + 2 * box, middle + box, SOFT, filled=False)
                start += 18 * s
            end = right
            if clear:
                self._text("see-through", right, middle, 10, SOFT, right=True)
                end -= 56 * s
            self._text(self._fit(label, end - start, 12), start, middle, 12, INK if on else SOFT)
            if toggle is not None:
                self.targets.append((left - 4 * s, middle - 10 * s, left + 17 * s, middle + 10 * s, toggle))
            self.targets.append((left, middle - 10 * s, right, middle + 10 * s, action))
            down(10)
            if kind is not None and session.palette == kind:
                palette()

        def palette():
            down(4)
            gap, size = _GAP * s, (right - left - 7 * _GAP * s) / 8
            for i, color in enumerate(PALETTE):
                x0 = left + (i % 8) * (size + gap)
                top = y[0] - (i // 8) * (size + gap)
                self._rect(x0, top - size, x0 + size, top, color)
                self.targets.append((x0, top - size, x0 + size, top,
                                     lambda color=color: session.pick_color(color)))
            y[0] -= 2 * size + 2 * gap
            buttons([("Default color", False, lambda: session.pick_color(None), 2)], 4)
            down(6)

        def item(kind, key, label, on, color, clear):
            check(label, on, lambda: session.row_clicked(kind, key), color, clear, (kind, key),
                  lambda: session.row_clicked(kind, key, "show"))

        if session.presets["title"]:
            self._text(self._fit(session.presets["title"], right - left, 17), left, down(10), 17, bold=True)
            down(14)
        for line in self._wrap(session.presets["note"], right - left, 11) if session.presets["note"] else ():
            self._text(line, left, down(7), 11, SOFT)
            down(8)
        heading("Cut along a section" if session.presets["sections"] else "View")
        current = session.current_section
        rows = [("Whole model", current is None, lambda: session.show_section(None), 4)]
        for row in session.presets["sections"]:
            rows.append((row["label"], current == row["label"],
                         lambda row=row: session.show_section(row["label"]),
                         2 if len(row["label"]) > 6 else 1))
        buttons(rows, 4)
        heading("Cut anywhere")
        buttons([(axis, session.slider_axis == axis and current == "", lambda axis=axis: session.cut_axis(axis), 1)
                 for axis in "xyz"] + [("Other side", False, session.flip_section, 2)], 5)
        self._slider(left, right, down(22))
        down(8)
        renderer = session.renderer
        heading("Colors and surfaces")
        buttons([("Custom colors", renderer.color_mode == "custom",
                  lambda: session.set_color_mode("custom"), 2),
                 ("Density colors", renderer.color_mode == "density",
                  lambda: session.set_color_mode("density"), 2),
                 ("All solid", renderer.surface == "solid", lambda: session.set_surfaces("solid"), 2),
                 ("All see-through", renderer.surface == "see-through",
                  lambda: session.set_surfaces("see-through"), 2)], 4)
        self._text("A click on a row below:", left, down(16), 11, SOFT)
        down(8)
        buttons([(label, session.tool == tool, lambda tool=tool: session.set_tool(tool), 1)
                 for label, tool in (("Show/hide", "show"), ("See-through", "clear"),
                                     ("Recolor", "color"))], 3)
        if session.anything_hidden():
            down(_GAP)
            buttons([("Show everything hidden", False, session.show_everything, 3)], 3)
        heading("Volumes")
        volumes = session.volume_rows()
        if not volumes:
            for line in self._wrap("Right-click a volume to list it here and change it on its own; "
                                   "shift+right-click hides it.",
                                   right - left, 11):
                self._text(line, left, down(7), 11, SOFT)
                down(8)
        for name, label, color, clear, on in volumes:
            item("volume", name, label, on, color, clear)
        if session.presets["layers"]:
            heading("Layers")
            names = [layer["label"] for layer in session.presets["layers"]]
            if OTHER in session.used_layers():
                names.append(OTHER)
            for name in names:
                item("layer", name, name, renderer.layers.get(name, True), *session.layer_row(name))
        materials = session.material_rows()
        if materials and session.legend:
            heading("Materials")
            for label, color, on, clear in materials:
                item("material", label, clean_name(label), on, color, clear)
        down(10)
        buttons([("Save settings", False, session.save_settings, 4)], 4)
        for line in self._wrap(session.save_message, right - left, 10) if session.save_message else ():
            self._text(line, left, down(7), 10, SOFT)
            down(7)
        self.height = (window_height + self.scroll - y[0]) + 14 * s

    def _tick(self, left, middle, box):
        points, cells = self.vtk.vtkPoints(), self.vtk.vtkCellArray()
        for x, y in ((.45, .05), (.85, -.4), (1.6, .5)):
            points.InsertNextPoint(left + x * box, middle + y * box, 0)
        cells.InsertNextCell(3)
        for i in range(3):
            cells.InsertCellPoint(i)
        data = self.vtk.vtkPolyData()
        data.SetPoints(points)
        data.SetLines(cells)
        mapper = self.vtk.vtkPolyDataMapper2D()
        mapper.SetInputData(data)
        actor = self.vtk.vtkActor2D()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1, 1, 1)
        actor.GetProperty().SetLineWidth(2 * self.scale)
        self._add(actor)

    def _slider(self, left, right, middle):
        s, session = self.scale, self.session
        end = right - 58 * s
        low, high = session.slider_range()
        value = min(max(session.slider_value, low), high)
        x = left + (end - left) * ((value - low) / (high - low) if high > low else .5)
        self._rect(left, middle - 1.5 * s, end, middle + 1.5 * s, RULE)
        self._rect(left, middle - 1.5 * s, x, middle + 1.5 * s, ACCENT)
        self._rect(x - 4 * s, middle - 8 * s, x + 4 * s, middle + 8 * s, ACCENT)
        self._text(("%.2f m" if abs(value) < 1e4 else "%.4g m") % value, right, middle, 11, SOFT, right=True)

        def slide(px, done=False):
            t = min(max((px - left) / (end - left), 0.), 1.)
            session.slide(low + t * (high - low), done)
        slide.draggable = True
        self.targets.append((left - 6 * s, middle - 11 * s, end + 6 * s, middle + 11 * s, slide))

    # ---- events ---------------------------------------------------------
    def press(self, x, y):
        """True when the panel handled a left press."""
        if not self.contains(x):
            return False
        for x0, y0, x1, y1, action in self.targets:
            if x0 <= x <= x1 and y0 <= y <= y1:
                if getattr(action, "draggable", False):
                    self.drag = action
                    action(x)
                else:
                    action()
                break
        return True

    def move(self, x):
        if self.drag is None:
            return False
        self.drag(x)
        return True

    def release(self, x):
        if self.drag is None:
            return False
        action, self.drag = self.drag, None
        action(x, True)
        return True

    def wheel(self, x, direction):
        if not self.contains(x):
            return False
        limit = max(self.height - self.window.GetSize()[1], 0.)
        self.scroll = min(max(self.scroll - direction * 40 * self.scale, 0.), limit)
        self.rebuild()
        return True
