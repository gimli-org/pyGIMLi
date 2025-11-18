#!/usr/bin/env python3
"""Not maintained anymore. Volunteers welcome to take over."""

import matplotlib as mpl
import numpy as np

# from PyQt5.QtCore import Qt, QPointF, QRect, QSize
# from PyQt5.QtGui import (
#     QPixmap, QPainter, QLinearGradient, QColor, QBrush, QDoubleValidator, QIcon
# )
# from PyQt5.QtWidgets import (
#     QWidget, QPushButton, QLineEdit, QComboBox, QSlider, QDoubleSpinBox,
#     QVBoxLayout, QHBoxLayout, QGroupBox, QLabel, QCheckBox
# )

QWidget = object
QPushButton = object
QLineEdit = object
QComboBox = object
QSlider = object
QDoubleSpinBox = object

# predefined color maps
CMAPS = ['viridis',
         'plasma',
         'Greys',
         'Blues',
         'Greens',
         'Reds',
         'gray',
         'hot',
         'Spectral',
         'seismic',
         'Set3',
         'terrain',
         'rainbow',
         'jet'
         ]


class GToolBar(QWidget):
    """
    Provide the toolbar for the 3D viewer with all buttons.
    This just provides the graphical features.
    """

    def __init__(self, parent=None):
        super(GToolBar, self).__init__(parent)
        self.setupWidget()

    def setupWidget(self):
        """
        Get some actions..
        """
        # combobox to choose the given parameter from
        self.cbbx_params = GComboBox("Scroll/Click to select parameter")

        # combobox to choose the colormap
        self.cbbx_cmap = GComboBox("Scroll/Click to select color map")
        for icon, name in self._createPixmap():
            self.cbbx_cmap.addItem(icon, name)
            self.cbbx_cmap.setIconSize(QSize(40, 15))

        # checkable button to reverse the chosen color scheme
        self.btn_reverse = GButton(
            text="Reverse",
            tooltip="Reverse the chosen color map",
            checkable=True,
            # size=[24, 24]
        )

        # # combobox for choosing the view
        # self.cbbx_views = GComboBox("Scroll/Click to select plot")
        # # checkbox for toggle link_views
        # self.chbx_views = QCheckBox("link views")

        # button to make bounding box visible
        self.btn_bbox = GButton(
            text="Toggle Axes",
            tooltip="Toggle data axis grid",
            checkable=True
        )
        self.btn_bbox.setChecked(True)

        # cMin and cMax
        self.spbx_cmin = GDoubleSpinBox("The min of the current range")
        self.spbx_cmin.setEnabled(False)
        self.spbx_cmax = GDoubleSpinBox("The max of the current range")
        self.spbx_cmax.setEnabled(False)

        # if not hit by key his button accepts changes to the color range
        self.btn_apply = GButton(
            text="Apply",
            tooltip="Apply changes in color range"
        )
        self.btn_apply.setEnabled(False)

        # resets the color range to its original range
        self.btn_reset = GButton(
            text="Reset",
            tooltip="Reset changes in color range"
        )
        self.btn_reset.setEnabled(False)

        # slider for slicing
        self.slice_x = GSlider("The X location of the YZ slice")
        self.slice_x.setEnabled(False)
        self.slice_y = GSlider("The Y location of the XZ slice")
        self.slice_y.setEnabled(False)
        self.slice_z = GSlider("The Z location of the XY slice")
        self.slice_z.setEnabled(False)

        # button to take a screenshot
        self.btn_screenshot = GButton(
            text="Screenshot",
            tooltip="Save screenshot of the current scene"
        )

        self.btn_exportVTK = GButton(
            text="Export VTK",
            tooltip="Save displayed mesh as VTK"
        )

        # parameter choosing
        lyt_params = QVBoxLayout()
        lyt_params.addWidget(self.cbbx_params)
        lyt_params.setContentsMargins(2, 5, 2, 2)
        grp_param = QGroupBox("Parameter")
        grp_param.setLayout(lyt_params)

        # lyt_views = QVBoxLayout()
        # lyt_views.setContentsMargins(2, 5, 2, 2)
        # lyt_views.addWidget(self.cbbx_views)
        # lyt_views.addWidget(self.chbx_views)
        # grp_views = QGroupBox("View")
        # grp_views.setLayout(lyt_views)

        # colormap
        lyt_cmap = QVBoxLayout()
        lyt_cmap.addWidget(self.cbbx_cmap)
        lyt_cmap.addWidget(self.btn_reverse)
        lyt_cmap.setContentsMargins(2, 5, 2, 2)
        grp_cmap = QGroupBox("Color Map")
        grp_cmap.setLayout(lyt_cmap)

        # limits
        lyt_min = QHBoxLayout()
        lyt_min.addWidget(QLabel("Min"))
        lyt_min.addWidget(self.spbx_cmin)
        lyt_max = QHBoxLayout()
        lyt_max.addWidget(QLabel("Max"))
        lyt_max.addWidget(self.spbx_cmax)
        # button for global limits
        self.btn_global_limits = GButton(
            text="Global Limits",
            tooltip="Check if gobal limits should be tried",
            checkable=True)
        # checkbox for thresholding parameter distribution
        lyt_thresh = QHBoxLayout()
        lyt_thresh.addWidget(QLabel("Threshold"))
        self.chbx_threshold = QCheckBox()
        lyt_thresh.addWidget(self.chbx_threshold)
        lyt_threshbtns = QHBoxLayout()
        lyt_threshbtns.addWidget(self.btn_apply)
        lyt_threshbtns.addWidget(self.btn_reset)
        lyt_threshold = QVBoxLayout()
        lyt_threshold.addLayout(lyt_min)
        lyt_threshold.addLayout(lyt_max)
        lyt_threshold.addWidget(self.btn_global_limits)
        lyt_threshold.addLayout(lyt_thresh)
        lyt_threshold.addLayout(lyt_threshbtns)
        lyt_threshold.setContentsMargins(2, 5, 2, 2)
        grp_limits = QGroupBox("Limits")
        grp_limits.setLayout(lyt_threshold)

        # slicer
        lyt_slicer = QVBoxLayout()
        self.btn_slice_plane = GButton(text="Plane", checkable=True)
        self.btn_slice_volume = GButton(text="Volume", checkable=True)
        # for the time being...
        self.btn_slice_volume.setEnabled(False)
        hb = QHBoxLayout()
        hb.addWidget(self.btn_slice_plane)
        hb.addWidget(self.btn_slice_volume)
        hx = QHBoxLayout()
        hx.addWidget(QLabel("x:"))
        hx.addWidget(self.slice_x)
        self.la_xval = QLineEdit()
        self.la_xval.setReadOnly(True)
        self.la_xval.setFixedWidth(100)
        hx.addWidget(self.la_xval)
        hy = QHBoxLayout()
        hy.addWidget(QLabel("y:"))
        hy.addWidget(self.slice_y)
        self.la_yval = QLineEdit()
        self.la_yval.setReadOnly(True)
        self.la_yval.setFixedWidth(100)
        hy.addWidget(self.la_yval)
        hz = QHBoxLayout()
        hz.addWidget(QLabel("z:"))
        hz.addWidget(self.slice_z)
        self.la_zval = QLineEdit()
        self.la_zval.setReadOnly(True)
        self.la_zval.setFixedWidth(100)
        hz.addWidget(self.la_zval)
        lyt_slicer.addLayout(hb)
        lyt_slicer.addLayout(hx)
        lyt_slicer.addLayout(hy)
        lyt_slicer.addLayout(hz)
        lyt_slicer.setContentsMargins(2, 5, 2, 2)
        self.grp_slice = QGroupBox("Slicing")
        # self.grp_slice.setCheckable(True)
        # self.grp_slice.setChecked(False)
        self.grp_slice.setLayout(lyt_slicer)

        # export area
        lyt_export = QHBoxLayout()
        lyt_export.addWidget(self.btn_screenshot)
        lyt_export.addWidget(self.btn_exportVTK)
        lyt_export.setContentsMargins(2, 5, 2, 2)
        grp_export = QGroupBox("Export")
        grp_export.setLayout(lyt_export)

        # widget for the toolbar to better organize the widgets
        lt = QVBoxLayout()
        lt.addWidget(grp_param)
        # lt.addWidget(grp_views)
        lt.addWidget(grp_limits)
        lt.addWidget(grp_cmap)
        lt.addWidget(self.grp_slice)
        lt.addStretch(1)
        lt.addWidget(self.btn_bbox)
        lt.addWidget(grp_export)
        lt.setContentsMargins(0, 0, 0, 0)

        self.setLayout(lt)

    def addExtraCMap(self, cMap):
        for icon, name in self._createPixmap([cMap]):
            self.cbbx_cmap.addItem(icon, name)
            self.cbbx_cmap.setIconSize(QSize(40, 15))

    def _createPixmap(self, cmaps=CMAPS):
        """
        Generator to create icons shown in the comboboxes where the colormap
        is selected.
        """
        w = 60
        h = 20
        for cmap in cmaps:
            cMap = mpl.cm.get_cmap(cmap)

            px = QPixmap(w, h)
            p = QPainter(px)
            gradient = QLinearGradient(QPointF(0, 0), QPointF(w, 0))

            for t in np.linspace(0., 1., 15):
                rgba = cMap(t)
                gradient.setColorAt(
                    t, QColor(rgba[0]*255, rgba[1]*255, rgba[2]*255))

            brush = QBrush(gradient)
            p.fillRect(QRect(0, 0, w, h), brush)
            p.end()

            yield QIcon(px), cmap


class GButton(QPushButton):

    def __init__(self, text=None, tooltip=None, checkable=False, size=None):
        super(GButton, self).__init__(None)
        self.setText(text)
        self.setToolTip(tooltip)
        self.setCheckable(checkable)
        if size is not None:
            self.setFixedSize(*size)


class GLineEdit(QLineEdit):

    def __init__(self, tooltip=None):
        super(GLineEdit, self).__init__(None)
        self.setToolTip(tooltip)
        self.setFixedWidth(200)
        # restrict acceptance to numbers only in that range
        # NOTE: will accept numbers outside of that range but the the
        # 'editingFinished' signal (hit enter/return when finished) won't
        self.setValidator(QDoubleValidator(-1e9, 1e9, 2, self))


class GComboBox(QComboBox):

    def __init__(self, tooltip=None):
        super(GComboBox, self).__init__(None)
        self.setToolTip(tooltip)


class _GDoubleSlider(QSlider):
    # https://gist.github.com/dennis-tra/994a65d6165a328d4eabaadbaedac2cc

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.decimals = 5
        self._max_int = 10 ** self.decimals

        super().setMinimum(0)
        super().setMaximum(self._max_int)

        self._min_value = -1e99
        self._max_value = 1e99

    @property
    def _value_range(self):
        return self._max_value - self._min_value

    def value(self):
        return float(super().value()) / self._max_int * self._value_range + self._min_value

    def setValue(self, value):
        super().setValue(int((value - self._min_value) / self._value_range * self._max_int))

    def setMinimum(self, value):
        if value > self._max_value:
            raise ValueError("Minimum limit cannot be higher than maximum")

        self._min_value = value
        self.setValue(self.value())

    def setMaximum(self, value):
        if value < self._min_value:
            raise ValueError("Minimum limit cannot be higher than maximum")

        self._max_value = value
        self.setValue(self.value())

    def minimum(self):
        return self._min_value

    def maximum(self):
        return self._max_value


class GSlider(_GDoubleSlider):

    def __init__(self, tooltip=None):
        super(GSlider, self).__init__(None)
        self.setOrientation(Qt.Horizontal)
        self.setToolTip(tooltip)


class GDoubleSpinBox(QDoubleSpinBox):

    def __init__(self, tooltip, parent=None):
        super(GDoubleSpinBox, self).__init__(parent)
        self.setToolTip(tooltip)
        self.setFixedWidth(200)
        self.setDecimals(6)
        self.setSingleStep(0.005)
        # self.setRange(*drange)
        self.setDecimals(10)
