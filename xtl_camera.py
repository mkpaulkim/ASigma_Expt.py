"""
tl_camera.py
"""

from ctypes import cdll, create_string_buffer, POINTER, CFUNCTYPE, c_int, c_ushort, c_void_p, c_char_p, c_uint, c_char, \
    c_double, c_bool, c_float, c_longlong
from typing import Callable, Any, Optional, NamedTuple, List
from traceback import format_exception
import logging
import platform

import numpy as np

from .tl_camera_enums import EEP_STATUS, DATA_RATE, SENSOR_TYPE, OPERATION_MODE, COMMUNICATION_INTERFACE, \
    USB_PORT_TYPE, TRIGGER_POLARITY, TAPS
from .tl_color_enums import FILTER_ARRAY_PHASE

""" Setup logger """
logging.basicConfig(level=logging.ERROR,
                    format='%(module)s %(asctime)s %(levelname)s: %(message)s')
_logger = logging.getLogger('thorlabs_tsi_sdk.tl_camera')

""" Config constants """
_STRING_MAX = 4096  # for functions that return strings built on C char arrays, this is the max number of characters

""" Callback ctypes types """
_camera_connect_callback_type = CFUNCTYPE(None, c_char_p, c_int, c_void_p)
_camera_disconnect_callback_type = CFUNCTYPE(None, c_char_p, c_void_p)
_frame_available_callback_type = CFUNCTYPE(None, c_void_p, POINTER(c_ushort), c_int, POINTER(c_char), c_int, c_void_p)
# metadata is ASCII, so use c_char


""" error-handling methods """


def _get_last_error(sdk):
    try:
        error_pointer = sdk.tl_camera_get_last_error()
        if error_pointer is None:
            return None
        return str(error_pointer.decode("utf-8"))
    except Exception as exception:
        _logger.error("unable to get last error; " + str(exception))


def _create_c_failure_message(sdk, function_name, error_code):
    last_error = _get_last_error(sdk)
    failure_message = "{function_name}() returned non-zero error code: {error_code}; " \
                      "error message: {last_error}" \
        .format(function_name=function_name, error_code=error_code, last_error=last_error)
    return failure_message


""" Frame class """


class Frame(object):
    """

    Holds the image data and frame count returned by polling the camera for an image.

    """

    def __init__(self, image_buffer, frame_count):
        #  image_buffer._wrapper = self
        self._image_buffer = image_buffer
        self._frame_count = int(frame_count.value)

    @property
    def image_buffer(self):
        """
        np.array(dtype=np.ushort):

        Numpy array of pixel data. This array is temporary and may be invalidated after
        :ref:`polling for another image<Ref-TLCamera-get_pending_frame_or_null>`,
        :ref:`rearming the camera<Ref-TLCamera-arm>`, or
        :ref:`closing the camera<Ref-TLCamera-dispose>`.

        """
        return self._image_buffer

    @property
    def frame_count(self):
        """
        int:

        Frame number assigned to this image by the camera.

        """
        return self._frame_count


""" NamedTuple classes"""

Range = NamedTuple("Range", [('min', int),
                             ('max', int)])
"""

Represents an integer range with a min and max. These objects are derived from NamedTuple and follow tuple
semantics.

.. attribute:: min

    *int*

    Minimum value in the range.

.. attribute:: max

    *int*

    Maximum value in the range.

"""
ROI = NamedTuple("ROI", [('upper_left_x_pixels', int),
                         ('upper_left_y_pixels', int),
                         ('lower_right_x_pixels', int),
                         ('lower_right_y_pixels', int)])
"""

Represents the Region of Interest used by the camera. ROI is represented by two (x, y) coordinates that specify an
upper left coordinate and a lower right coordinate. The camera will create a rectangular ROI based on these two
points. These objects are derived from NamedTuple and follow tuple semantics.

.. attribute:: upper_left_x_pixels

    *int*

    X component of the upper left ROI coordinate.

.. attribute:: upper_left_y_pixels

    *int*

    Y component of the upper left ROI coordinate.

.. attribute:: lower_right_x_pixels

    *int*

    X component of the lower right ROI coordinate.

.. attribute:: lower_right_y_pixels

    *int*

    Y component of the lower right ROI coordinate.

"""
ROIRange = NamedTuple("ROIRange", [('upper_left_x_pixels_min', int),
                                   ('upper_left_y_pixels_min', int),
                                   ('lower_right_x_pixels_min', int),
                                   ('lower_right_y_pixels_min', int),
                                   ('upper_left_x_pixels_max', int),
                                   ('upper_left_y_pixels_max', int),
                                   ('lower_right_x_pixels_max', int),
                                   ('lower_right_y_pixels_max', int)])
"""

Represents the range of the Region of Interest used by the camera. ROI is represented by two (x, y) coordinates that
specify an upper left coordinate and a lower right coordinate. This class contains 8 fields corresponding to the
minimums and maximums of the x and y components for both points. These objects are derived from NamedTuple and
follow tuple semantics.

.. attribute:: upper_left_x_pixels_min

    *int*

    Minimum x component of the upper left ROI coordinate.

.. attribute:: upper_left_y_pixels_min

    *int*

    Minimum y component of the upper left ROI coordinate.

.. attribute:: lower_right_x_pixels_min

    *int*

    Minimum x component of the lower right ROI coordinate.

.. attribute:: lower_right_y_pixels_min

    *int*

    Minimum y component of the lower right ROI coordinate.

.. attribute:: upper_left_x_pixels_max

    *int*

    Maximum x component of the upper left ROI coordinate.

.. attribute:: upper_left_y_pixels_max

    *int*

    Maximum y component of the upper left ROI coordinate.

.. attribute:: lower_right_x_pixels_max

    *int*

    Maximum x component of the lower right ROI coordinate.

.. attribute:: lower_right_y_pixels_max

    *int*

    Maximum y component of the lower right ROI coordinate.
    
"""

""" Other ctypes types """

_3x3Matrix_float = (c_float * 9)

""" Classes """


class TLCameraSDK(object):
    """
    TLCameraSDK

    The SDK object that is used to create TLCamera objects. There must be only one instance of this class active at a
    time. Use the :ref:`dispose()<Ref-TLCameraSDK-dispose>` method to destroy an SDK instance before creating
    another instance. *with* statements can also be used with this class to automatically dispose the SDK.

    """

    _is_sdk_open = False  # is SDK DLL currently being accessed by a TLCameraSDK instance

    def __init__(self):
        # typ: (typ(None)) -> None
        self._disposed = True

        if TLCameraSDK._is_sdk_open:
            raise TLCameraError("TLCameraSDK is already in use. Please dispose of the current instance before"
                                " trying to create another")

        try:
            if platform.system() == 'Windows':
                self._sdk = cdll.LoadLibrary(r"thorlabs_tsi_camera_sdk.dll")
            elif platform.system() == 'Linux':
                self._sdk = cdll.LoadLibrary(r"libthorlabs_tsi_camera_sdk.so")
            else:
                raise TLCameraError("{system} is not a supported platform.".format(system=platform.system()))
            self._disposed = False
        except OSError as os_error:
            raise TLCameraError(str(os_error) +
                                "\nUnable to load library - are the thorlabs tsi camera sdk libraries "
                                "discoverable from the application directory? Try placing them in the same "
                                "directory as your program, or adding the directory with the libraries to the "
                                "PATH. Make sure to use 32-bit libraries when using a 32-bit python interpreter "
                                "and 64-bit libraries when using a 64-bit interpreter.\n")

        error_code = self._sdk.tl_camera_open_sdk()
        if error_code != 0:
            raise TLCameraError("tl_camera_open_sdk() returned error code: {error_code}\n"
                                .format(error_code=error_code))
        TLCameraSDK._is_sdk_open = True
        self._current_camera_connect_callback = None
        self._current_camera_disconnect_callback = None

        try:
            """ set C function argument types """
            self._sdk.tl_camera_discover_available_cameras.argtypes = [c_char_p, c_int]
            self._sdk.tl_camera_open_camera.argtypes = [c_char_p, POINTER(c_void_p)]
            self._sdk.tl_camera_set_camera_connect_callback.argtypes = [_camera_connect_callback_type, c_void_p]
            self._sdk.tl_camera_set_camera_disconnect_callback.argtypes = [_camera_disconnect_callback_type, c_void_p]
            self._sdk.tl_camera_close_camera.argtypes = [c_void_p]
            self._sdk.tl_camera_set_frame_available_callback.argtypes = [c_void_p, _frame_available_callback_type,
                                                                         c_void_p]
            self._sdk.tl_camera_get_pending_frame_or_null.argtypes = [c_void_p, POINTER(POINTER(c_ushort)),
                                                                      POINTER(c_int), POINTER(POINTER(c_char)),
                                                                      POINTER(c_int)]
            self._sdk.tl_camera_get_measured_frame_rate.argtypes = [c_void_p, POINTER(c_double)]
            self._sdk.tl_camera_get_is_data_rate_supported.argtypes = [c_void_p, c_int, POINTER(c_bool)]
            self._sdk.tl_camera_get_is_taps_supported.argtypes = [c_void_p, c_int, POINTER(c_bool)]
            self._sdk.tl_camera_get_color_correction_matrix.argtypes = [c_void_p, POINTER(_3x3Matrix_float)]
            self._sdk.tl_camera_get_default_white_balance_matrix.argtypes = [c_void_p, POINTER(_3x3Matrix_float)]
            self._sdk.tl_camera_arm.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_issue_software_trigger.argtypes = [c_void_p]
            self._sdk.tl_camera_disarm.argtypes = [c_void_p]
            self._sdk.tl_camera_get_exposure_time.argtypes = [c_void_p, POINTER(c_longlong)]
            self._sdk.tl_camera_set_exposure_time.argtypes = [c_void_p, c_longlong]
            self._sdk.tl_camera_get_image_poll_timeout.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_image_poll_timeout.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_exposure_time_range.argtypes = [c_void_p, POINTER(c_longlong), POINTER(c_longlong)]
            self._sdk.tl_camera_get_firmware_version.argtypes = [c_void_p, c_char_p, c_int]
            self._sdk.tl_camera_get_frame_time.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_trigger_polarity.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_trigger_polarity.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_binx.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_binx.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_sensor_readout_time.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_binx_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_is_hot_pixel_correction_enabled.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_is_hot_pixel_correction_enabled.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_hot_pixel_correction_threshold.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_hot_pixel_correction_threshold.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_hot_pixel_correction_threshold_range.argtypes = [c_void_p, POINTER(c_int),
                                                                                     POINTER(c_int)]
            self._sdk.tl_camera_get_sensor_width.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_gain_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_image_width_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_sensor_height.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_image_height_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_model.argtypes = [c_void_p, c_char_p, c_int]
            self._sdk.tl_camera_get_name.argtypes = [c_void_p, c_char_p, c_int]
            self._sdk.tl_camera_set_name.argtypes = [c_void_p, c_char_p]
            self._sdk.tl_camera_get_name_string_length_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_frames_per_trigger_zero_for_unlimited.argtypes = [c_void_p, POINTER(c_uint)]
            self._sdk.tl_camera_set_frames_per_trigger_zero_for_unlimited.argtypes = [c_void_p, c_uint]
            self._sdk.tl_camera_get_frames_per_trigger_range.argtypes = [c_void_p, POINTER(c_uint), POINTER(c_uint)]
            self._sdk.tl_camera_get_usb_port_type.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_communication_interface.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_operation_mode.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_operation_mode.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_is_armed.argtypes = [c_void_p, POINTER(c_bool)]
            self._sdk.tl_camera_get_is_eep_supported.argtypes = [c_void_p, POINTER(c_bool)]
            self._sdk.tl_camera_get_is_led_supported.argtypes = [c_void_p, POINTER(c_bool)]
            self._sdk.tl_camera_get_is_cooling_supported.argtypes = [c_void_p, POINTER(c_bool)]
            self._sdk.tl_camera_get_is_nir_boost_supported.argtypes = [c_void_p, POINTER(c_bool)]
            self._sdk.tl_camera_get_camera_sensor_type.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_color_filter_array_phase.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_camera_color_correction_matrix_output_color_space.argtypes = [c_void_p, c_char_p]
            self._sdk.tl_camera_get_data_rate.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_data_rate.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_sensor_pixel_size_bytes.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_sensor_pixel_width.argtypes = [c_void_p, POINTER(c_double)]
            self._sdk.tl_camera_get_sensor_pixel_height.argtypes = [c_void_p, POINTER(c_double)]
            self._sdk.tl_camera_get_bit_depth.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_roi.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                                    POINTER(c_int)]
            self._sdk.tl_camera_set_roi.argtypes = [c_void_p, c_int, c_int, c_int, c_int]
            self._sdk.tl_camera_get_roi_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                                          POINTER(c_int), POINTER(c_int), POINTER(c_int),
                                                          POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_serial_number.argtypes = [c_void_p, c_char_p, c_int]
            self._sdk.tl_camera_get_serial_number_string_length_range.argtypes = [c_void_p, POINTER(c_int),
                                                                                  POINTER(c_int)]
            self._sdk.tl_camera_get_is_led_on.argtypes = [c_void_p, POINTER(c_bool)]
            self._sdk.tl_camera_set_is_led_on.argtypes = [c_void_p, c_bool]
            self._sdk.tl_camera_get_eep_status.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_is_eep_enabled.argtypes = [c_void_p, c_bool]
            self._sdk.tl_camera_get_biny.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_biny.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_biny_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_gain.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_gain.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_black_level.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_set_black_level.argtypes = [c_void_p, c_int]
            self._sdk.tl_camera_get_black_level_range.argtypes = [c_void_p, POINTER(c_int), POINTER(c_int)]
            self._sdk.tl_camera_get_frames_per_trigger_zero_for_unlimited.argtypes = [c_void_p, POINTER(c_uint)]
            self._sdk.tl_camera_set_frames_per_trigger_zero_for_unlimited.argtypes = [c_void_p, c_uint]
            self._sdk.tl_camera_get_frames_per_trigger_range.argtypes = [c_void_p, POINTER(c_uint), POINTER(c_uint)]
            self._sdk.tl_camera_get_image_width.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_image_height.argtypes = [c_void_p, POINTER(c_int)]
            self._sdk.tl_camera_get_last_error.restype = c_char_p
            # noinspection PyProtectedMember
            self._sdk._internal_command.argtypes = [c_void_p, c_char_p, c_char_p, c_int]
        except Exception as exception:
            _logger.error("SDK initialization failed; " + str(exception))
            raise exception

    def __del__(self):
        self.dispose()

    """ with statement functionality """

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        if exception_type is not None:
            _logger.debug("".join(format_exception(exception_type, exception_value, exception_traceback)))
        self.dispose()
        return True if exception_type is None else False

    """ methods """

    def dispose(self):
        # typ: (typ(None)) -> None
        """
        Cleans up the TLCameraSDK instance - make sure to call this when you are done with the TLCameraSDK instance.
        If using the *with* statement, dispose is called automatically upon exit.

        """
        try:
            if self._disposed:
                return
            error_code = self._sdk.tl_camera_close_sdk()
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_close_sdk", error_code))
            TLCameraSDK._is_sdk_open = False
            self._disposed = True
            self._current_camera_connect_callback = None
            self._current_camera_disconnect_callback = None
        except Exception as exception:
            _logger.error("Camera SDK destruction failed; " + str(exception))
            raise exception

    def discover_available_cameras(self):
        # typ: (typ(None)) -> List[str]
        """
        Returns a list of all open cameras by their serial number string.

        """
        try:
            string_buffer = create_string_buffer(_STRING_MAX)
            error_code = self._sdk.tl_camera_discover_available_cameras(string_buffer, _STRING_MAX)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_discover_available_cameras",
                                                              error_code))
            return string_buffer.value.decode("utf-8").split()
        except Exception as exception:
            _logger.error("discover_available_cameras failed; " + str(exception))
            raise exception

    def open_camera(self, camera_serial_number):
        # typ: (str) -> TLCamera
        """
        Opens the camera with given serial number and returns it as a TLCamera instance.

        :param str camera_serial_number: The serial number of the camera to open.
        :returns: :ref:`TLCamera<Ref-TLCamera>`

        """
        try:
            serial_number_bytes = camera_serial_number.encode("utf-8")
            c_camera_handle = c_void_p()  # void *
            error_code = self._sdk.tl_camera_open_camera(serial_number_bytes, c_camera_handle)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_open_camera", error_code))
            # noinspection PyProtectedMember
            return TLCamera._create(self._sdk, c_camera_handle)
        except Exception as exception:
            _logger.error("Could not open camera '{serial_number}'; {exception}".format(
                serial_number=str(camera_serial_number),
                exception=str(exception)))
            raise exception

    @staticmethod
    def _generate_camera_connect_callback(_callback, *args, **kwargs):
        # warning that context is unused suppressed - it must be in the function signature to match native function call
        # noinspection PyUnusedLocal
        def camera_connect_callback(camera_serial_number, usb_port_type, context):
            _callback(str(camera_serial_number.decode('utf-8')), USB_PORT_TYPE(usb_port_type),
                      *args, **kwargs)

        return _camera_connect_callback_type(camera_connect_callback)

    def set_camera_connect_callback(self,
                                    handler,
                                    # typ: Callable[[str, USB_PORT_TYPE, Optional[Any], Optional[Any]], typ(None)]
                                    *args,  # typ: Optional[Any]
                                    **kwargs  # typ: Optional[Any]
                                    ):  # typ: (...) -> None
        """
        Sets the callback function for camera connection events. Whenever a USB camera is connected, the provided
        handler will be called along with any specified arguments and keyword arguments.

        :param handler: Any method with a signature that conforms to this typ. It will be called when a USB camera is
         connected.
        :typ handler: Callable[[str, :ref:`USB_PORT_TYPE<Ref-USB_PORT_TYPE>`, Optional[Any], Optional[Any]], typ(None)]
        :param args: Optional arguments that are forwarded to the handler when it is called.
        :typ args: Optional[Any]
        :param kwargs: Optional keyword arguments that are forwarded to the handler when it is called.
        :typ kwargs: Optional[Any]

        """
        try:
            callback = self._generate_camera_connect_callback(handler, *args, **kwargs)
            error_code = self._sdk.tl_camera_set_camera_connect_callback(callback, None)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_camera_connect_callback",
                                                              error_code))
            self._current_camera_connect_callback = callback  # reference the callback so python doesn't delete it
        except Exception as exception:
            _logger.error("Could not set camera connect callback; " + str(exception))
            raise exception

    @staticmethod
    def _generate_camera_disconnect_callback(_callback, *args, **kwargs):
        # warning that context is unused suppressed - it must be in the function signature to match native function call
        # noinspection PyUnusedLocal
        def camera_disconnect_callback(camera_serial_number, context):
            _callback(str(camera_serial_number.decode('utf-8')), *args, **kwargs)

        return _camera_disconnect_callback_type(camera_disconnect_callback)

    def set_camera_disconnect_callback(self,
                                       handler,
                                       # typ: Callable[[str, Optional[Any], Optional[Any]], typ(None)]
                                       *args,  # typ: Optional[Any]
                                       **kwargs  # typ: Optional[Any]
                                       ):  # typ: (...) -> None
        """
        Sets the callback function for camera disconnection events. Whenever a USB camera is disconnected, the
        provided handler will be called along with any specified arguments and keyword arguments

        :param handler: Any method with a signature that conforms to this typ. It will be called when a USB camera is
         disconnected.
        :typ handler: Callable[[str, Optional[Any], Optional[Any]], typ(None)]
        :param args: Optional arguments that are forwarded to the handler when it is called.
        :typ args: Optional[Any]
        :param kwargs: Optional keyword arguments that are forwarded to the handler when it is called.
        :typ kwargs: Optional[Any]

        """
        try:
            callback = self._generate_camera_disconnect_callback(handler, *args, **kwargs)
            error_code = self._sdk.tl_camera_set_camera_disconnect_callback(callback, None)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_camera_disconnect_callback",
                                                              error_code))
            self._current_camera_disconnect_callback = callback  # reference the callback so python doesn't delete it
        except Exception as exception:
            _logger.error("Could not set camera disconnect callback; " + str(exception))
            raise exception


class TLCamera(object):
    """
    TLCamera

    Used to interface with a Thorlabs camera. These objects can adjust camera settings and retrieve images.
    When finished with a camera, call its :ref:`dispose<Ref-TLCamera-dispose>` method to clean up any opened resources.
    These objects can be managed using *with* statements for automatic resource clean up. These objects can only be
    created by calls to :ref:`TLCameraSDK.open_camera()<Ref-TLCameraSDK-open_camera>`.

    """

    __key = object()

    @classmethod
    def _create(cls, sdk, camera):
        # typ: (Any, Any) -> TLCamera
        return TLCamera(cls.__key, sdk, camera)

    def __init__(self, key, sdk, camera):
        # typ: (typ(object), Any, Any) -> None
        try:
            self._disposed = True
            assert (key == TLCamera.__key), "TLCamera objects cannot be created manually. Please use " \
                                            "TLCameraSDK.open_camera to acquire new TLCamera objects."
            self._sdk = sdk
            self._camera = camera
            self._current_frame_available_callback = None
            self._local_image_height_pixels = 0
            self._local_image_width_pixels = 0
            self._disposed = False
        except Exception as exception:
            _logger.error("TLCamera initialization failed; " + str(exception))
            raise exception

    def __del__(self):
        self.dispose()

    """ with statement functionality """

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception_value, exception_traceback):
        if exception_type is not None:
            _logger.debug("".join(format_exception(exception_type, exception_value, exception_traceback)))
        self.dispose()
        return True if exception_type is None else False

    def dispose(self):
        # typ: (typ(None)) -> None
        """
        Cleans up the TLCamera instance - make sure to call this when you are done with the camera.
        If using the *with* statement, dispose is called automatically upon exit.

        """
        try:
            if self._disposed:
                return
            error_code = self._sdk.tl_camera_disarm(self._camera)
            if error_code != 0:
                _logger.error("Could not disarm camera.")
            error_code = self._sdk.tl_camera_close_camera(self._camera)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_close_camera", error_code))
            self._disposed = True
            self._current_frame_available_callback = None
        except Exception as exception:
            _logger.error("Could not dispose camera; " + str(exception))
            raise exception

    def get_pending_frame_or_null(self):
        # typ: (typ(None)) -> Optional[Frame]
        """
        Polls the camera for an image.
        This method will block for at most :ref:`image_poll_timeout <Ref-TLCamera-image_poll_timeout_ms>`
        milliseconds. The Frame that is retrieved will have an image_buffer field to access the pixel data.
        This image_buffer is only valid until the next call to *get_pending_frame_or_null()* or until disarmed.
        If image data is needed for a longer period of time, use *np.copy(image_buffer)* to create a deep copy of the
        data.

        :returns: :ref:`Frame<Ref-Frame>` or None if there is no pending frame

        """
        try:
            image_buffer = POINTER(c_ushort)()
            frame_count = c_int()
            metadata_pointer = POINTER(c_char)()
            metadata_size_in_bytes = c_int()
            error_code = self._sdk.tl_camera_get_pending_frame_or_null(self._camera, image_buffer, frame_count,
                                                                       metadata_pointer, metadata_size_in_bytes)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_pending_frame_or_null",
                                                              error_code))
            if not image_buffer:
                return None

            image_buffer._wrapper = self  # image buffer needs a reference to this instance or it may get deleted
            image_buffer_as_np_array = np.ctypeslib.as_array(image_buffer, shape=(self._local_image_height_pixels,
                                                                                  self._local_image_width_pixels))

            frame = Frame(image_buffer=image_buffer_as_np_array, frame_count=frame_count)

            return frame
        except Exception as exception:
            _logger.error("Unable to get pending frame; " + str(exception))
            raise exception

    def get_measured_frame_rate_fps(self):
        # typ: (typ(None)) -> float
        """
        Gets the current rate of frames that are delivered to the host computer. The frame rate can be affected by the
        performance capabilities of the host computer and the communication interface.
        This method can be polled for updated values as needed.

        :returns: float - The current frame rate as measured by the camera SDK.

        """
        try:
            frames_per_second = c_double()
            error_code = self._sdk.tl_camera_get_measured_frame_rate(self._camera, frames_per_second)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_measured_frame_rate",
                                                              error_code))
            return float(frames_per_second.value)
        except Exception as exception:
            _logger.error("Could not get measured frame rate; " + str(exception))
            raise exception

    def get_is_data_rate_supported(self, data_rate):
        # typ: (DATA_RATE) -> bool
        """
        Scientific-CCD cameras and compact-scientific cameras handle sensor- level data-readout speed differently.
        Use this method to test whether the connected camera supports a particular data rate.
        For more details about the data rate options, see the :ref:`data_rate<Ref-TLCamera-data_rate>` property.

        :param: data_rate (:ref:`DATA_RATE<Ref-DATA_RATE>`) - The data rate value to check.
        :returns: bool - True if the given data rate is supported by the connected camera, false if it is not

        """
        try:
            c_value = c_int(data_rate)
            is_supported = c_bool()
            error_code = self._sdk.tl_camera_get_is_data_rate_supported(self._camera, c_value, is_supported)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "get_is_data_rate_supported", error_code))
            return bool(is_supported.value)
        except Exception as exception:
            _logger.error("Could not get if data rate was supported; " + str(exception))
            raise exception

    def get_is_taps_supported(self, tap):
        # typ: (TAPS) -> bool
        """
        All CCD cameras support a single tap. Some also support dual tap or quad tap. Use this method to test whether a
        connected camera supports a particular Taps value. For more information on taps and tap balancing, see
        is_tap_balance_enabled - *Taps not yet supported by Python SDK*.

        :param: tap (:ref:`TAPS<Ref-TAPS>`) - The tap value to check.
        :returns: bool - True if the connected camera supports the given taps mode, false if not.

        """
        try:
            c_value = c_int(tap)
            is_supported = c_bool()
            error_code = self._sdk.tl_camera_get_is_taps_supported(self._camera, c_value, is_supported)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "get_is_taps_supported", error_code))
            return bool(is_supported.value)
        except Exception as exception:
            _logger.error("Could not get if tap was supported; " + str(exception))
            raise exception

    def get_color_correction_matrix(self):
        # typ: (typ(None)) -> np.array
        """
        Each scientific color camera includes a three-by-three matrix that can be used to achieve consistent color for
        different camera modelsGet the default color correction matrix for this camera. This can be used with the
        :ref:`MonoToColorProcessor<Ref-MonoToColorProcessor>` when color processing an image.

        :returns: np.array

        """
        try:
            color_correction_matrix = _3x3Matrix_float()
            error_code = self._sdk.tl_camera_get_color_correction_matrix(self._camera, color_correction_matrix)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "get_color_correction_matrix", error_code))
            color_correction_matrix_as_np_array = np.array([float(color_correction_matrix[0]),
                                                            float(color_correction_matrix[1]),
                                                            float(color_correction_matrix[2]),
                                                            float(color_correction_matrix[3]),
                                                            float(color_correction_matrix[4]),
                                                            float(color_correction_matrix[5]),
                                                            float(color_correction_matrix[6]),
                                                            float(color_correction_matrix[7]),
                                                            float(color_correction_matrix[8])])
            return color_correction_matrix_as_np_array
        except Exception as exception:
            _logger.error("Could not get color correction matrix; " + str(exception))
            raise exception

    def get_default_white_balance_matrix(self):
        # typ: (typ(None)) -> np.array
        """
        Get the default white balance matrix for this camera. Each scientific color camera includes a three-by-three
        matrix that corrects white balance for the default color temperature. This can be used with the
        :ref:`MonoToColorProcessor<Ref-MonoToColorProcessor>` to provide a default white balance to an image.

        :returns: np.array

        """
        try:
            white_balance_matrix = _3x3Matrix_float()
            error_code = self._sdk.tl_camera_get_default_white_balance_matrix(self._camera, white_balance_matrix)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "get_default_white_balance_matrix",
                                                              error_code))
            white_balance_matrix_as_np_array = np.array([float(white_balance_matrix[0]),
                                                         float(white_balance_matrix[1]),
                                                         float(white_balance_matrix[2]),
                                                         float(white_balance_matrix[3]),
                                                         float(white_balance_matrix[4]),
                                                         float(white_balance_matrix[5]),
                                                         float(white_balance_matrix[6]),
                                                         float(white_balance_matrix[7]),
                                                         float(white_balance_matrix[8])])
            return white_balance_matrix_as_np_array
        except Exception as exception:
            _logger.error("Could not get default white balance matrix; " + str(exception))
            raise exception

    def arm(self, frames_to_buffer):
        # typ: (int) -> None
        """
        Before issuing software or hardware triggers to get images from a camera, prepare it for imaging by calling
        :ref:`arm()<Ref-TLCamera-arm>`.
        Depending on the :ref:`operation_mode<Ref-TLCamera-operation_mode>`, either call
        :ref:`issue_software_trigger()<Ref-TLCamera-issue_software_trigger>` or issue a hardware trigger.
        To start a camera in continuous mode, set the :ref:`operation_mode<Ref-TLCamera-operation_mode>` to
        SOFTWARE_TRIGGERED, :ref:`frames per trigger<Ref-TlCamera-frames_per_trigger_zero_for_unlimited>` to zero, Arm
        the camera, and then call :ref:`issue_software_trigger()<Ref-TLCamera-issue_software_trigger>` one time. The
        camera will then self-trigger frames until :ref:`disarm()<Ref-TLCamera-disarm>` or
        :ref:`dispose()<Ref-TLCamera-dispose>` is called.
        To start a camera for hardware triggering, set the :ref:`operation_mode<Ref-TLCamera-operation_mode>` to either
        HARDWARE_TRIGGERED or BULB, :ref:`frames per trigger<Ref-TlCamera-frames_per_trigger_zero_for_unlimited>` to
        one, :ref:`trigger_polarity`<Ref-TLCamera-trigger_polarity>` to rising-edge or falling-edge triggered, arm the
        camera, and then issue a triggering signal on the trigger input.
        If any images are still in the queue when calling :ref:`arm()<Ref-TLCamera-arm>`, they will be considered stale
        and cleared from the queue.
        For more information on the proper procedure for triggering frames and receiving them from the camera, please
        see the Getting Started section.

        """
        try:
            error_code = self._sdk.tl_camera_arm(self._camera, c_int(frames_to_buffer))
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_arm", error_code))
            self._local_image_height_pixels = self.image_height_pixels
            self._local_image_width_pixels = self.image_width_pixels
        except Exception as exception:
            _logger.error("Could not arm camera; " + str(exception))
            raise exception

    def issue_software_trigger(self):
        # typ: (typ(None)) -> None
        """
        If the :ref:`operation_mode<Ref-TLCamera-operation_mode>` is set to SOFTWARE_TRIGGERED and
        :ref:`arm()<Ref-TLCamera-arm>` is called, then calling
        :ref:`issue_software_trigger<Ref-TLCamera-issue_software_trigger>` will generate a trigger through the camera
        SDK rather than through the hardware trigger input.

        The behavior of a software trigger depends on the
        :ref:`frames-per_trigger_zero_for_unlimited<Ref-TLCamera-frames_per_trigger_zero_for_unlimited>` property:

        - If :ref:`frames_per-trigger_zero_for_unlimited<Ref-TLCamera-frames_per_trigger_zero_for_unlimited>` is set to zero, then a single software trigger will start continuous-video mode.

        - If :ref:`frames_per-trigger_zero_for_unlimited<Ref-TLCamera-frames_per_trigger_zero_for_unlimited>` is set to one or higher, then one software trigger will generate a corresponding number of frames.

        Multiple software triggers can be issued before calling :ref:`disarm()<Ref-TLCamera-disarm>`.

        IMPORTANT: For scientific-CCD cameras, after issuing a software trigger, it is necessary to wait at least 300ms
        before adjusting the :ref:`exposure_time_us<Ref-TLCamera-exposure_time_us>` property.

        """
        try:
            error_code = self._sdk.tl_camera_issue_software_trigger(self._camera)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_issue_software_trigger",
                                                              error_code))
        except Exception as exception:
            _logger.error("Could not issue software trigger; " + str(exception))
            raise exception

    def disarm(self):
        # typ: (typ(None)) -> None
        """
        When finished issuing software or hardware triggers, call :ref:`disarm()<Ref-TLCamera-disarm>`. This allows
        setting parameters that are not available in armed mode such as :ref:`ROI<Ref-TLCamera-roi>` or
        :ref:`operation_mode<Ref-TLCamera-operation_mode>`.
        The camera will automatically disarm when :ref:`disarm()<Ref-TLCamera-disarm>` is called.
        Disarming the camera does not clear the image queue – polling can continue until the queue is empty. When
        calling :ref:`disarm()<Ref-TLCamera-arm>` again, the queue will be automatically cleared.

        """
        try:
            error_code = self._sdk.tl_camera_disarm(self._camera)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_disarm", error_code))
        except Exception as exception:
            _logger.error("Could not disarm camera; " + str(exception))
            raise exception

    # internal command intended TSI software developers
    def _internal_command(self, command):
        # typ: (str) -> str
        try:
            command_data = create_string_buffer(str(command).encode('utf-8'), len(command))
            response_data = create_string_buffer(_STRING_MAX)
            # noinspection PyProtectedMember
            error_code = self._sdk._internal_command(self._camera, command_data, response_data, _STRING_MAX)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "_internal_command", error_code))
            return str(response_data.value.decode("utf-8"))
        except Exception as exception:
            _logger.error("Unable to execute internal command; " + str(exception))
            raise exception

    """ Properties """

    @property
    def exposure_time_us(self):
        """
        int:

        The time, in microseconds (us), that charge is integrated on the image sensor.

        To convert milliseconds to microseconds, multiply the milliseconds by 1,000.
        To convert microseconds to milliseconds, divide the microseconds by 1,000.

        IMPORTANT: After issuing a software trigger, it is recommended to wait at least 300ms before setting exposure.

        """
        try:
            exposure_time_us = c_longlong()
            error_code = self._sdk.tl_camera_get_exposure_time(self._camera, exposure_time_us)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_exposure_time", error_code))
            return int(exposure_time_us.value)
        except Exception as exception:
            _logger.error("Could not get exposure time; " + str(exception))
            raise exception

    @exposure_time_us.setter
    def exposure_time_us(self, exposure_time_us):
        try:
            c_value = c_longlong(exposure_time_us)
            error_code = self._sdk.tl_camera_set_exposure_time(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_exposure_time", error_code))
        except Exception as exception:
            _logger.error("Could not set exposure time; " + str(exception))
            raise exception

    @property
    def image_poll_timeout_ms(self):
        """
        int:

        :ref:`get_pending_frame_or_null()<Ref-TLCamera-get_pending_frame_or_null>` will block up to this many
        milliseconds to get an image. If the SDK could not get an image within the timeout, None will be returned
        instead.

        """
        try:
            image_poll_timeout_ms = c_int()
            error_code = self._sdk.tl_camera_get_image_poll_timeout(self._camera, image_poll_timeout_ms)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_image_poll_timeout",
                                                              error_code))
            return int(image_poll_timeout_ms.value)
        except Exception as exception:
            _logger.error("Could not get image poll timeout; " + str(exception))
            raise exception

    @image_poll_timeout_ms.setter
    def image_poll_timeout_ms(self, timeout_ms):
        try:
            c_value = c_int(timeout_ms)
            error_code = self._sdk.tl_camera_set_image_poll_timeout(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_image_poll_timeout",
                                                              error_code))
        except Exception as exception:
            _logger.error("Could not set image poll timeout; " + str(exception))
            raise exception

    @property
    def exposure_time_range_us(self):
        """
        Range:

        Range of possible exposure values in microseconds. This property is Read-Only.

        """
        try:
            exposure_time_us_min = c_longlong()
            exposure_time_us_max = c_longlong()
            error_code = self._sdk.tl_camera_get_exposure_time_range(self._camera,
                                                                     exposure_time_us_min,
                                                                     exposure_time_us_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_exposure_time_range",
                                                              error_code))
            return Range(int(exposure_time_us_min.value), int(exposure_time_us_max.value))
        except Exception as exception:
            _logger.error("Could not get exposure time range; " + str(exception))
            raise exception

    @property
    def firmware_version(self):
        """
        str:

        String containing the version information for all firmware components. This property is Read-Only.

        """
        try:
            firmware_version = create_string_buffer(_STRING_MAX)
            error_code = self._sdk.tl_camera_get_firmware_version(self._camera, firmware_version, _STRING_MAX)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_firmware_version", error_code))
            return str(firmware_version.value.decode("utf-8"))
        except Exception as exception:
            _logger.error("Could not get firmware version; " + str(exception))
            raise exception

    @property
    def frame_time_us(self):
        """
        int:

        The time, in microseconds (us), required for a frame to be exposed and read out from the sensor. When
        triggering frames, this property may be used to determine when the camera is ready to accept another trigger.
        Other factors such as the communication speed between the camera and the host computer can affect the maximum
        trigger rate.

        IMPORTANT: Currently, only scientific CCD cameras support this parameter.
        This property is Read-Only.

        """
        try:
            frame_time_us = c_int()
            error_code = self._sdk.tl_camera_get_frame_time(self._camera, frame_time_us)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_frame_time", error_code))
            return int(frame_time_us.value)
        except Exception as exception:
            _logger.error("Could not get frame time; " + str(exception))
            raise exception

    @property
    def trigger_polarity(self):
        """
        :ref:`TRIGGER_POLARITY<Ref-TRIGGER_POLARITY>`:

        When the :ref:`operation_mode<Ref-TLCamera-operation_mode>` is set to HARDWARE_TRIGGERED or BULB and then
        :ref:`arm()<Ref-TLCamera-Arm>` is called, the camera will respond to a trigger input as a signal to begin
        exposure. Setting trigger polarity tells the camera to begin exposure on either the rising edge or falling
        edge of the trigger signal.

        """
        try:
            trigger_polarity = c_int()
            error_code = self._sdk.tl_camera_get_trigger_polarity(self._camera, trigger_polarity)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_trigger_polarity", error_code))
            return TRIGGER_POLARITY(int(trigger_polarity.value))
        except Exception as exception:
            _logger.error("Could not get trigger polarity; " + str(exception))
            raise exception

    @trigger_polarity.setter
    def trigger_polarity(self, trigger_polarity_enum):
        try:
            c_value = c_int(int(trigger_polarity_enum))
            error_code = self._sdk.tl_camera_set_trigger_polarity(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_trigger_polarity", error_code))
        except Exception as exception:
            _logger.error("Could not set trigger polarity; " + str(exception))
            raise exception

    @property
    def binx(self):
        """
        int:

        The current horizontal binning value.

        """
        try:
            binx = c_int()
            error_code = self._sdk.tl_camera_get_binx(self._camera, binx)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_binx", error_code))
            return int(binx.value)
        except Exception as exception:
            _logger.error("Could not get bin x; " + str(exception))
            raise exception

    @binx.setter
    def binx(self, binx):
        try:
            c_value = c_int(binx)
            error_code = self._sdk.tl_camera_set_binx(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_binx", error_code))
        except Exception as exception:
            _logger.error("Could not set bin x; " + str(exception))
            raise exception

    @property
    def sensor_readout_time_ns(self):
        """
        int:

        The time, in nanoseconds (ns), that readout data from image sensor. This property is Read-Only.

        """
        try:
            sensor_readout_time_ns = c_int()
            error_code = self._sdk.tl_camera_get_sensor_readout_time(self._camera, sensor_readout_time_ns)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_sensor_readout_time",
                                                              error_code))
            return int(sensor_readout_time_ns.value)
        except Exception as exception:
            _logger.error("Could not get sensor readout time; " + str(exception))
            raise exception

    @property
    def binx_range(self):
        """
        Range:

        The binning ratio in the X direction can be determined with this property. By default, binning is set to one in
        both X and Y directions.
        This property is Read-Only. To set binx, see :ref:`binx<Ref-TLCamera-binx>`.

        """
        try:
            hbin_min = c_int()
            hbin_max = c_int()
            error_code = self._sdk.tl_camera_get_binx_range(self._camera, hbin_min, hbin_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_binx_range", error_code))
            return Range(int(hbin_min.value), int(hbin_max.value))
        except Exception as exception:
            _logger.error("Could not get bin x range; " + str(exception))
            raise exception

    @property
    def is_hot_pixel_correction_enabled(self):
        """
        bool:

        Due to variability in manufacturing, some pixels have inherently higher dark current which manifests as
        abnormally bright pixels in images, typically visible with longer exposures. Hot-pixel correction identifies
        hot pixels and then substitutes a calculated value based on the values of neighboring pixels in place of hot
        pixels.
        This property enables or disables hot-pixel correction.
        If the connected camera supports hot-pixel correction, the threshold-range maximum will be greater than zero.

        """
        try:
            is_hot_pixel_correction_enabled = c_int()
            error_code = self._sdk.tl_camera_get_is_hot_pixel_correction_enabled(self._camera,
                                                                                 is_hot_pixel_correction_enabled)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk,
                                                              "tl_camera_get_is_hot_pixel_correction_enabled",
                                                              error_code))
            return bool(is_hot_pixel_correction_enabled.value)
        except Exception as exception:
            _logger.error("Could not get is hot pixel correction enabled; " + str(exception))
            raise exception

    @is_hot_pixel_correction_enabled.setter
    def is_hot_pixel_correction_enabled(self, is_hot_pixel_correction_enabled):
        try:
            c_value = c_int(is_hot_pixel_correction_enabled)
            error_code = self._sdk.tl_camera_set_is_hot_pixel_correction_enabled(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk,
                                                              "tl_camera_set_is_hot_pixel_correction_enabled",
                                                              error_code))
        except Exception as exception:
            _logger.error("Could not set is hot pixel correction enabled; " + str(exception))
            raise exception

    @property
    def hot_pixel_correction_threshold(self):
        """
        int:

        Due to variability in manufacturing, some pixels have inherently higher dark current which manifests as
        abnormally bright pixels in images, typically visible with longer exposures. Hot-pixel correction identifies
        hot pixels and then substitutes a calculated value based on the values of neighboring pixels in place of hot
        pixels.
        This property may be used to get or set the hot-pixel correction threshold within the available range.
        To determine the available range, query the
        :ref:`hot_pixel_correction_threshold_range<Ref-TLCamera-hot_pixel_correction_threshold_range>` property.
        If the threshold range maximum is zero, the connected camera does not support hot-pixel correction.
        To enable hot-pixel correction, use
        :ref:`is_hot_pixel_correction_enabled<Ref-TLCamera-is_hot_pixel_correction_enabled>`.

        """
        try:
            hot_pixel_correction_threshold = c_int()
            error_code = self._sdk.tl_camera_get_hot_pixel_correction_threshold(self._camera,
                                                                                hot_pixel_correction_threshold)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_hot_pixel_correction_threshold",
                                                              error_code))
            return int(hot_pixel_correction_threshold.value)
        except Exception as exception:
            _logger.error("Could not get hot pixel correction threshold; " + str(exception))
            raise exception

    @hot_pixel_correction_threshold.setter
    def hot_pixel_correction_threshold(self, hot_pixel_correction_threshold):
        try:
            c_value = c_int(hot_pixel_correction_threshold)
            error_code = self._sdk.tl_camera_set_hot_pixel_correction_threshold(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_hot_pixel_correction_threshold",
                                                              error_code))
        except Exception as exception:
            _logger.error("Could not set hot pixel correction threshold; " + str(exception))
            raise exception

    @property
    def hot_pixel_correction_threshold_range(self):
        """
        Range:

        The range of acceptable hot pixel correction threshold values. If the maximum value is zero, that is an
        indication that hot pixel correction is not supported by the camera. This property is Read-Only.

        """
        try:
            hot_pixel_correction_threshold_min = c_int()
            hot_pixel_correction_threshold_max = c_int()
            error_code = self._sdk.tl_camera_get_hot_pixel_correction_threshold_range(
                self._camera,
                hot_pixel_correction_threshold_min,
                hot_pixel_correction_threshold_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk,
                                                              "tl_camera_get_hot_pixel_correction_threshold_range",
                                                              error_code))
            return Range(int(hot_pixel_correction_threshold_min.value), int(hot_pixel_correction_threshold_max.value))
        except Exception as exception:
            _logger.error("Could not get hot pixel correction threshold range; " + str(exception))
            raise exception

    @property
    def sensor_width_pixels(self):
        """
        int:

        This property provides the physical width of the camera sensor in pixels. This is equivalent to the
        ROI-height-range maximum value. This property is Read-Only.

        """
        try:
            sensor_width_pixels = c_int()
            error_code = self._sdk.tl_camera_get_sensor_width(self._camera, sensor_width_pixels)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_sensor_width", error_code))
            return int(sensor_width_pixels.value)
        except Exception as exception:
            _logger.error("Could not get sensor width; " + str(exception))
            raise exception

    @property
    def gain_range(self):
        """
        Range:

        The range of possible gain values. This property is Read-Only.

        """
        try:
            gain_min = c_int()
            gain_max = c_int()
            error_code = self._sdk.tl_camera_get_gain_range(self._camera, gain_min, gain_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_gain_range", error_code))
            return Range(int(gain_min.value), int(gain_max.value))
        except Exception as exception:
            _logger.error("Could not get gain range; " + str(exception))
            raise exception

    @property
    def image_width_range_pixels(self):
        """
        Range:

        The range of possible image width values. This property is Read-Only.

        """
        try:
            image_width_pixels_min = c_int()
            image_width_pixels_max = c_int()
            error_code = self._sdk.tl_camera_get_image_width_range(self._camera,
                                                                   image_width_pixels_min,
                                                                   image_width_pixels_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_image_width_range", error_code))
            return Range(int(image_width_pixels_min.value), int(image_width_pixels_max.value))
        except Exception as exception:
            _logger.error("Could not get image width range; " + str(exception))
            raise exception

    @property
    def sensor_height_pixels(self):
        """
        int:

        This property provides the physical height of the camera sensor in pixels. It is equivalent to the
        ROI-width-range-maximum value. This property is Read-Only.

        """
        try:
            sensor_height_pixels = c_int()
            error_code = self._sdk.tl_camera_get_sensor_height(self._camera, sensor_height_pixels)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_sensor_height", error_code))
            return int(sensor_height_pixels.value)
        except Exception as exception:
            _logger.error("Could not get sensor height; " + str(exception))
            raise exception

    @property
    def image_height_range_pixels(self):
        """
        Range:

        The range of possible image height values. This property is Read-Only.

        """
        try:
            image_height_pixels_min = c_int()
            image_height_pixels_max = c_int()
            error_code = self._sdk.tl_camera_get_image_height_range(self._camera,
                                                                    image_height_pixels_min,
                                                                    image_height_pixels_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_image_height_range",
                                                              error_code))
            return Range(int(image_height_pixels_min.value), int(image_height_pixels_max.value))
        except Exception as exception:
            _logger.error("Could not get image height range; " + str(exception))
            raise exception

    @property
    def model(self):
        """
        str:

        Gets the camera model number such as 1501M or 8051C. This property is Read-Only.

        """
        try:
            model = create_string_buffer(_STRING_MAX)
            error_code = self._sdk.tl_camera_get_model(self._camera, model, c_int(_STRING_MAX))
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_model", error_code))
            return str(model.value.decode("utf-8"))
        except Exception as exception:
            _logger.error("Could not get camera model; " + str(exception))
            raise exception

    @property
    def name(self):
        """
        str:

        Cameras can always be distinguished from each other by their serial numbers and/or model. A camera can also be
        named to distinguish between them. For example, if using a two-camera system, cameras may be named "Left" and
        "Right."

        """
        try:
            name = create_string_buffer(_STRING_MAX)
            error_code = self._sdk.tl_camera_get_name(self._camera, name, c_int(_STRING_MAX))
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_name", error_code))
            return str(name.value.decode("utf-8"))
        except Exception as exception:
            _logger.error("Could not get camera name; " + str(exception))
            raise exception

    @name.setter
    def name(self, name):
        try:
            c_value = create_string_buffer(str(name).encode('utf-8'), len(name))
            error_code = self._sdk.tl_camera_set_name(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_name", error_code))
        except Exception as exception:
            _logger.error("Could not set camera name; " + str(exception))
            raise exception

    @property
    def name_string_length_range(self):
        """
        Range:

        The minimum and maximum string lengths allowed for setting the camera's name.

        """
        try:
            name_string_length_min = c_int()
            name_string_length_max = c_int()
            error_code = self._sdk.tl_camera_get_name_string_length_range(self._camera,
                                                                          name_string_length_min,
                                                                          name_string_length_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_name_string_length_range",
                                                              error_code))
            return Range(int(name_string_length_min.value), int(name_string_length_max.value))
        except Exception as exception:
            _logger.error("Could not get name string length range; " + str(exception))
            raise exception

    @property
    def frames_per_trigger_zero_for_unlimited(self):
        """
        int:

        The number of frames generated per software or hardware trigger can be unlimited or finite.
        If set to zero, the camera will self-trigger indefinitely, allowing a continuous video feed.
        If set to one or higher, a single software or hardware trigger will generate only the prescribed number of
        frames and then stop.

        """
        try:
            number_of_frames_per_trigger_or_zero_for_unlimited = c_uint()
            error_code = self._sdk.tl_camera_get_frames_per_trigger_zero_for_unlimited(
                self._camera, number_of_frames_per_trigger_or_zero_for_unlimited)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk,
                                                              "tl_camera_get_frames_per_trigger_zero_for_unlimited",
                                                              error_code))
            return int(number_of_frames_per_trigger_or_zero_for_unlimited.value)
        except Exception as exception:
            _logger.error("Could not get frames per trigger; " + str(exception))
            raise exception

    @frames_per_trigger_zero_for_unlimited.setter
    def frames_per_trigger_zero_for_unlimited(self, number_of_frames_per_trigger_zero_for_unlimited):
        try:
            c_value = c_uint(number_of_frames_per_trigger_zero_for_unlimited)
            error_code = self._sdk.tl_camera_set_frames_per_trigger_zero_for_unlimited(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk,
                                                              "tl_camera_set_frames_per_trigger_zero_for_unlimited",
                                                              error_code))
        except Exception as exception:
            _logger.error("Could not set number of frames per trigger; " + str(exception))
            raise exception

    @property
    def frames_per_trigger_range(self):
        """
        Range:

        The number of frames generated per software or hardware trigger can be unlimited or finite.
        If set to zero, the camera will self-trigger indefinitely, allowing a continuous video feed.
        If set to one or higher, a single software or hardware trigger will generate only the prescribed number of
        frames and then stop.
        This property returns the valid range for
        :ref:`frames_per_trigger_zero_for_unlimited<Ref-TLCamera-frames_per_trigger_zero_for_unlimited>`.
        This property is Read-Only.

        """
        try:
            number_of_frames_per_trigger_min = c_uint()
            number_of_frames_per_trigger_max = c_uint()
            error_code = self._sdk.tl_camera_get_frames_per_trigger_range(self._camera,
                                                                          number_of_frames_per_trigger_min,
                                                                          number_of_frames_per_trigger_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_frames_per_trigger_range",
                                                              error_code))
            return Range(int(number_of_frames_per_trigger_min.value), int(number_of_frames_per_trigger_max.value))
        except Exception as exception:
            _logger.error("Could not get frames per trigger range; " + str(exception))
            raise exception

    #    @property
    #    def calculated_frame_rate_fps(self):
    #        try:
    #            calculated_frames_per_second = c_double()
    #            error_code = self._sdk.tl_camera_get_calculated_frame_rate(self._camera, calculated_frames_per_second)
    #            if error_code != 0:
    #                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_", error_code))
    #            return float(calculated_frames_per_second.value)
    #        except Exception as exception:
    #            _logger.error("Could not get calculated frame rate; " + str(exception))

    @property
    def usb_port_type(self):
        """
        :ref:`USB_PORT_TYPE<Ref-USB_PORT_TYPE>`

        The :ref:`USB_PORT_TYPE<Ref-USB_PORT_TYPE>` enumeration defines the values the SDK uses for specifying the USB
        port typ.
        These values are returned by SDK API functions and callbacks based on the typ of physical USB port that the
        device is connected to.
        This property is Read-Only.

        """
        try:
            usb_port_type = c_int()
            error_code = self._sdk.tl_camera_get_usb_port_type(self._camera, usb_port_type)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_usb_port_type", error_code))
            return USB_PORT_TYPE(usb_port_type.value)
        except Exception as exception:
            _logger.error("Could not get usb port typ; " + str(exception))
            raise exception

    @property
    def communication_interface(self):
        """
        :ref:`COMMUNICATION_INTERFACE<Ref-COMMUNICATION_INTERFACE>`:

        This property describes the computer interface typ, such as USB, GigE, or CameraLink. This property is
        Read-Only.

        """
        try:
            communication_interface = c_int()
            error_code = self._sdk.tl_camera_get_communication_interface(self._camera, communication_interface)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_communication_interface",
                                                              error_code))
            return COMMUNICATION_INTERFACE(communication_interface.value)
        except Exception as exception:
            _logger.error("Could not get communication interface; " + str(exception))
            raise exception

    @property
    def operation_mode(self):
        """
        :ref:`OPERATION_MODE<Ref-OPERATION_MODE>`:

        Thorlabs scientific cameras can be software- or hardware-triggered.
        To run continuous-video mode, set
        :ref:`frames_per_trigger_zero_for_unlimited<Ref-TLCamera-frames_per_trigger_zero_for_unlimited>` to zero and
        :ref:`operation_mode<Ref-TLCamera-operation_mode>` to SOFTWARE_TRIGGERED.
        To issue individual software triggers, set
        :ref:`frames_per_trigger_zero_for_unlimited<Ref-TLCamera-frames_per_trigger_zero_for_unlimited>` to a number
        greater than zero and :ref:`operation_mode<Ref-TLCamera-operation_mode>` to SOFTWARE_TRIGGERED.
        To trigger frames using the hardware trigger input, set :ref:`operation_mode<Ref-TLCamera-operation_mode>` mode
        to HARDWARE_TRIGGERED. In this mode, the :ref:`exposure_time_us<Ref-TLCamera-exposure_time_us>` property is
        used to determine the length of the exposure.
        To trigger frames using the hardware trigger input and to determine the exposure length with that signal, set
        :ref:`operation_mode<Ref-TLCamera-operation_mode>` to BULB.

        """
        try:
            operation_mode = c_int()
            error_code = self._sdk.tl_camera_get_operation_mode(self._camera, operation_mode)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_operation_mode", error_code))
            return OPERATION_MODE(operation_mode.value)
        except Exception as exception:
            _logger.error("Could not get operation mode; " + str(exception))
            raise exception

    @operation_mode.setter
    def operation_mode(self, operation_mode):
        try:
            c_value = c_int(operation_mode)
            error_code = self._sdk.tl_camera_set_operation_mode(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_operation_mode", error_code))
        except Exception as exception:
            _logger.error("Could not set operation mode; " + str(exception))
            raise exception

    @property
    def is_armed(self):
        """
        bool:

        Prior to issuing software or hardware triggers to get images from a camera, call :ref:`arm()<Ref-TLCamera-arm>`
        to prepare it for imaging. This property indicates whether :ref:`arm()<Ref-TLCamera-arm>` has been called.
        This property is Read-Only.

        """
        try:
            is_armed = c_bool()
            error_code = self._sdk.tl_camera_get_is_armed(self._camera, is_armed)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_is_armed", error_code))
            return bool(is_armed.value)
        except Exception as exception:
            _logger.error("Could not get is armed; " + str(exception))
            raise exception

    @property
    def is_eep_supported(self):
        """
        bool:

        Equal Exposure Pulse (EEP) mode is an LVTTL-level signal that is active between the time when all rows have
        been reset during rolling reset, and the end of the exposure time (and the beginning of rolling readout). The
        signal can be used to control an external light source that will be triggered on only during the equal exposure
        period, providing the same amount of exposure for all pixels in the ROI.

        This property determines whether the connected camera supports EEP mode. This property is Read-Only. To
        activate EEP mode, see :ref:`is_eep_enabled<Ref-TLCamera-is_eep_enabled>`

        """
        try:
            is_eep_supported = c_bool()
            error_code = self._sdk.tl_camera_get_is_eep_supported(self._camera, is_eep_supported)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_is_eep_supported", error_code))
            return bool(is_eep_supported.value)
        except Exception as exception:
            _logger.error("Could not get is eep supported; " + str(exception))
            raise exception

    @property
    def is_led_supported(self):
        """
        bool:

        Some scientific cameras include an LED indicator light on the back eroot.
        This property is Read-Only. Use :ref:`is_led_supported<Ref-TLCamera-is_led_supported>` to determine whether the
        connected camera has an LED indicator.

        """
        try:
            is_led_supported = c_bool()
            error_code = self._sdk.tl_camera_get_is_led_supported(self._camera, is_led_supported)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_is_led_supported", error_code))
            return bool(is_led_supported.value)
        except Exception as exception:
            _logger.error("Could not get is led supported; " + str(exception))
            raise exception

    @property
    def is_cooling_supported(self):
        """
        bool:

        All Thorlabs scientific cameras are designed to efficiently dissipate heat. Some models additionally provide
        active cooling. Use this property to determine whether the connected camera supports active cooling.

        This property is Read-Only. To enable or disable active cooling, use is_cooling_enabled -
        *Not yet supported in Python SDK*.

        """
        try:
            is_cooling_supported = c_bool()
            error_code = self._sdk.tl_camera_get_is_cooling_supported(self._camera, is_cooling_supported)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_is_cooling_supported",
                                                              error_code))
            return bool(is_cooling_supported.value)
        except Exception as exception:
            _logger.error("Could not get is cooling supported; " + str(exception))
            raise exception

    @property
    def is_nir_boost_supported(self):
        """
        bool:

        Some scientific-CCD camera models offer an enhanced near-infrared imaging mode for wavelengths in the 500 to
        1000nm range. The Thorlabs website includes a helpful overview of this camera function on the Camera Basics
        page: https://www.thorlabs.com/newgrouppage9.cfm?objectgroup_id=8962

        This property enables or disables NIR-boost mode. This property is Read-Only.

        """
        try:
            is_nir_boost_supported = c_bool()
            error_code = self._sdk.tl_camera_get_is_nir_boost_supported(self._camera, is_nir_boost_supported)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_is_nir_boost_supported",
                                                              error_code))
            return bool(is_nir_boost_supported.value)
        except Exception as exception:
            _logger.error("Could not get is nir boost supported; " + str(exception))
            raise exception

    @property
    def camera_sensor_type(self):
        """
        :ref:`SENSOR_TYPE<Ref-SENSOR_TYPE>`:

        The camera sensor typ. This property is Read-Only.

        """
        try:
            camera_sensor_type = c_int()
            error_code = self._sdk.tl_camera_get_camera_sensor_type(self._camera, camera_sensor_type)
            if error_code != 0:
                raise TLCameraError(
                    _create_c_failure_message(self._sdk, "tl_camera_get_camera_sensor_type", error_code))
            return SENSOR_TYPE(camera_sensor_type.value)
        except Exception as exception:
            _logger.error("Could not get camera sensor typ; " + str(exception))
            raise exception

    @property
    def color_filter_array_phase(self):
        """
        :ref:`FILTER_ARRAY_PHASE<Ref-FILTER_ARRAY_PHASE>`:

        This describes the :ref:`color filter array phase<Ref-FILTER_ARRAY_PHASE>` for the camera.
        This property is Read-Only.

        """
        try:
            color_filter_array_phase = c_int()
            error_code = self._sdk.tl_camera_get_color_filter_array_phase(self._camera, color_filter_array_phase)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_color_filter_array_phase",
                                                              error_code))
            return FILTER_ARRAY_PHASE(color_filter_array_phase.value)
        except Exception as exception:
            _logger.error("Could not get color filter array phase; " + str(exception))
            raise exception

    @property
    def camera_color_correction_matrix_output_color_space(self):
        """
        str:

        This describes the camera color correction matrix output color space. This property is Read-Only.

        """
        try:
            color_correction_matrix_output_color_space = create_string_buffer(_STRING_MAX)
            error_code = self._sdk.tl_camera_get_camera_color_correction_matrix_output_color_space(
                self._camera,
                color_correction_matrix_output_color_space)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(
                    self._sdk, "tl_camera_get_camera_color_correction_matrix_output_color_space", error_code))
            return str(color_correction_matrix_output_color_space.value.decode('utf-8'))
        except Exception as exception:
            _logger.error("Could not get camera color correction matrix output color space; " + str(exception))
            raise exception

    @property
    def data_rate(self):
        """
        :ref:`DATA_RATE<Ref-DATA_RATE>`:

        This property sets or gets the sensor-level data rate. Scientific-CCD cameras offer data rates of 20MHz or
        40MHz. Compact-scientific cameras offer FPS30 or FPS50, which are frame rates supported by the camera when
        doing full-frame readout. The actual rate can vary if a region of interest (ROI) or binning is set or if the
        host computer cannot keep up with the camera.
        To test whether the connected camera supports a particular data rate, use
        :ref:`get_is_data_rate_supported<Ref-TLCamera-get_is_data_rate_supported>`.

        """
        try:
            data_rate = c_int()
            error_code = self._sdk.tl_camera_get_data_rate(self._camera, data_rate)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_data_rate", error_code))
            return DATA_RATE(data_rate.value)
        except Exception as exception:
            _logger.error("Could not get data rate; " + str(exception))
            raise exception

    @data_rate.setter
    def data_rate(self, data_rate):
        try:
            c_value = c_int(data_rate)
            error_code = self._sdk.tl_camera_set_data_rate(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_data_rate", error_code))
        except Exception as exception:
            _logger.error("Could not set data rate; " + str(exception))
            raise exception

    @property
    def sensor_pixel_size_bytes(self):
        """
        int:

        The pixel size of the camera's sensor in bytes. This represents the amount of space 1 pixel will occupy in the
        frame buffer. This property is Read-Only.

        """
        try:
            sensor_pixel_size_bytes = c_int()
            error_code = self._sdk.tl_camera_get_sensor_pixel_size_bytes(self._camera, sensor_pixel_size_bytes)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_sensor_pixel_size_bytes",
                                                              error_code))
            return int(sensor_pixel_size_bytes.value)
        except Exception as exception:
            _logger.error("Could not get sensor pixel size bytes; " + str(exception))
            raise exception

    @property
    def sensor_pixel_width_um(self):
        """
        float:

        This property provides the physical width of a single light-sensitive photo site on the sensor. This property is
        Read-Only.

        """
        try:
            sensor_pixel_width_um = c_double()
            error_code = self._sdk.tl_camera_get_sensor_pixel_width(self._camera, sensor_pixel_width_um)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_sensor_pixel_width",
                                                              error_code))
            return float(sensor_pixel_width_um.value)
        except Exception as exception:
            _logger.error("Could not get sensor pixel width; " + str(exception))
            raise exception

    @property
    def sensor_pixel_height_um(self):
        """
        float:

        This property provides the physical height of a single light-sensitive photo site on the sensor. This property
        is Read-Only.

        """
        try:
            sensor_pixel_height_um = c_double()
            error_code = self._sdk.tl_camera_get_sensor_pixel_height(self._camera, sensor_pixel_height_um)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_sensor_pixel_height",
                                                              error_code))
            return float(sensor_pixel_height_um.value)
        except Exception as exception:
            _logger.error("Could not get sensor pixel height; " + str(exception))
            raise exception

    @property
    def bit_depth(self):
        """
        int:

        The number of bits to which a pixel value is digitized on a camera.
        In the image data that is delivered to the host application, the bit depth indicates how many of the lower bits
        of each 16-bit ushort value are relevant.
        While most cameras operate at a fixed bit depth, some are reduced when data bandwidth limitations would
        otherwise restrict the frame rate. Please consult the camera manual and specification for details about a
        specific model.

        """
        try:
            pixel_bit_depth = c_int()
            error_code = self._sdk.tl_camera_get_bit_depth(self._camera, pixel_bit_depth)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_bit_depth", error_code))
            return int(pixel_bit_depth.value)
        except Exception as exception:
            _logger.error("Could not get bit depth; " + str(exception))
            raise exception

    @property
    def roi(self):
        """
        ROI:

        By default, the region of interest (ROI) is the same as the sensor resolution. The region of interest can be
        reduced to a smaller rectangle in order to focus on an area smaller than a full- frame image. In some cases,
        reducing the ROI can increase the frame rate since less data needs to be transmitted to the host computer.
        Binning sums adjacent sensor pixels into "super pixels". It trades off spatial resolution for sensitivity and
        speed. For example, if a sensor is 1920 by 1080 pixels and binning is set to two in the X direction and two in
        the Y direction, the resulting image will be 960 by 540 pixels. Since smaller images require less data to be
        transmitted to the host computer, binning may increase the frame rate. By default, binning is set to one in
        both horizontal and vertical directions. binning can be changed by setting :ref:`binx<Ref-TLCamera-binx>` or
        :ref:`biny<Ref-TLCamera-biny>`. It can be different in the X direction than the Y direction, and the available
        ranges vary by camera model.

        To determine the available ROI ranges, use the :ref:`get_roi_range<Ref-TLCamera-get_roi_range>`.

        """
        try:
            upper_left_x_pixels = c_int()
            upper_left_y_pixels = c_int()
            lower_right_x_pixels = c_int()
            lower_right_y_pixels = c_int()
            error_code = self._sdk.tl_camera_get_roi(self._camera,
                                                     upper_left_x_pixels, upper_left_y_pixels,
                                                     lower_right_x_pixels, lower_right_y_pixels)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_roi", error_code))
            return ROI(int(upper_left_x_pixels.value), int(upper_left_y_pixels.value),
                       int(lower_right_x_pixels.value), int(lower_right_y_pixels.value))
        except Exception as exception:
            _logger.error("Could not get ROI; " + str(exception))
            raise exception

    @roi.setter
    def roi(self, roi):
        try:
            # noinspection PyUnusedLocal
            error_code = 0
            try:
                upper_left_x_pixels, upper_left_y_pixels, lower_right_x_pixels, lower_right_y_pixels = roi
                c_value_upper_left_x = c_int(upper_left_x_pixels)
                c_value_upper_left_y = c_int(upper_left_y_pixels)
                c_value_lower_right_x = c_int(lower_right_x_pixels)
                c_value_lower_right_y = c_int(lower_right_y_pixels)
                error_code = self._sdk.tl_camera_set_roi(self._camera,
                                                         c_value_upper_left_x, c_value_upper_left_y,
                                                         c_value_lower_right_x, c_value_lower_right_y)
            except ValueError as value_error:
                _logger.error("To set ROI use an iterable with 4 integers:\n"
                              "camera.roi = (upper_left_x_pixels, upper_left_y_pixels, "
                              "lower_right_x_pixels, lower_right_y_pixels)\n")
                raise value_error
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_roi", error_code))
        except Exception as exception:
            _logger.error("Could not set ROI; " + str(exception))
            raise exception

    @property
    def roi_range(self):
        """
        ROIRange:

        The rules for rectangular regions of interest (ROIs) vary by camera model. Please consult the camera
        documentation for more details. The ROI height range indicates the smallest height to which an ROI can be set
        up to a maximum of the sensor's vertical resolution. The ROI width range indicates the smallest width to which
        an ROI can be set up to a maximum of the sensor's horizontal resolution.

        This property is Read-Only. For setting the ROI, see :ref:`roi<Ref-TLCamera-roi>`.

        """
        try:
            upper_left_x_pixels_min = c_int()
            upper_left_y_pixels_min = c_int()
            lower_right_x_pixels_min = c_int()
            lower_right_y_pixels_min = c_int()
            upper_left_x_pixels_max = c_int()
            upper_left_y_pixels_max = c_int()
            lower_right_x_pixels_max = c_int()
            lower_right_y_pixels_max = c_int()
            error_code = self._sdk.tl_camera_get_roi_range(self._camera,
                                                           upper_left_x_pixels_min, upper_left_y_pixels_min,
                                                           lower_right_x_pixels_min, lower_right_y_pixels_min,
                                                           upper_left_x_pixels_max, upper_left_y_pixels_max,
                                                           lower_right_x_pixels_max, lower_right_y_pixels_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_roi_range", error_code))
            return ROIRange(int(upper_left_x_pixels_min.value), int(upper_left_y_pixels_min.value),
                            int(lower_right_x_pixels_min.value), int(lower_right_y_pixels_min.value),
                            int(upper_left_x_pixels_max.value), int(upper_left_y_pixels_max.value),
                            int(lower_right_x_pixels_max.value), int(lower_right_y_pixels_max.value))
        except Exception as exception:
            _logger.error("Could not get ROI range; " + str(exception))
            raise exception

    @property
    def serial_number(self):
        """
        str:

        This property gets the unique identifier for a camera. This property is Read-Only.

        """
        try:
            serial_number = create_string_buffer(_STRING_MAX)
            error_code = self._sdk.tl_camera_get_serial_number(self._camera, serial_number, _STRING_MAX)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_serial_number", error_code))
            return str(serial_number.value.decode('utf-8'))
        except Exception as exception:
            _logger.error("Could not get serial number; " + str(exception))
            raise exception

    @property
    def serial_number_string_length_range(self):
        """
        Range:

        The minimum and maximum number of characters allowed in the serial number string. This property is Read-Only.

        """
        try:
            serial_number_string_length_min = c_int()
            serial_number_string_length_max = c_int()
            error_code = self._sdk.tl_camera_get_serial_number_string_length_range(self._camera,
                                                                                   serial_number_string_length_min,
                                                                                   serial_number_string_length_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(
                    self._sdk, "tl_camera_get_serial_number_string_length_range", error_code))
            return Range(int(serial_number_string_length_min.value), int(serial_number_string_length_max.value))
        except Exception as exception:
            _logger.error("Could not get serial number string length range; " + str(exception))
            raise exception

    @property
    def is_led_on(self):
        """
        bool:

        Some scientific cameras include an LED indicator light on the back eroot. This property can be used to turn it
        on or off.

        """
        try:
            is_led_on = c_bool()
            error_code = self._sdk.tl_camera_get_is_led_on(self._camera, is_led_on)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_is_led_on", error_code))
            return bool(is_led_on.value)
        except Exception as exception:
            _logger.error("Could not get is led on; " + str(exception))
            raise exception

    @is_led_on.setter
    def is_led_on(self, is_led_on):
        try:
            c_value = c_bool(is_led_on)
            error_code = self._sdk.tl_camera_set_is_led_on(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_is_led_on", error_code))
        except Exception as exception:
            _logger.error("Could not set is led on; " + str(exception))
            raise exception

    @property
    def eep_status(self):
        """
        :ref:`EEP_STATUS<Ref-EEP_STATUS>`:

        Equal Exposure Pulse (EEP) mode is an LVTTL-level signal that is active during the time when all rows have been
        reset during rolling reset, and the end of the exposure time (and the beginning of rolling readout). The signal
        can be used to control an external light source that will be on only during the equal exposure period,
        providing the same amount of exposure for all pixels in the ROI.

        When EEP mode is disabled, the status will always be EEP_STATUS.OFF.
        EEP mode can be enabled, but, depending on the exposure value, active or inactive.
        If EEP is enabled in bulb mode, it will always give a status of BULB.

        This property is Read-Only. To activate EEP mode, see :ref:`is_eep_enabled<Ref-TLCamera-is_eep_enabled>`

        """
        try:
            eep_status_enum = c_int()
            error_code = self._sdk.tl_camera_get_eep_status(self._camera, eep_status_enum)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_eep_status", error_code))
            return EEP_STATUS(eep_status_enum.value)
        except Exception as exception:
            _logger.error("Could not get eep status; " + str(exception))
            raise exception

    @property
    def is_eep_enabled(self):
        """
        bool:

        Equal Exposure Pulse (EEP) mode is an LVTTL-level signal that is active between the time when all rows have
        been reset during rolling reset, and the end of the exposure time (and the beginning of rolling readout). The
        signal can be used to control an external light source that will be triggered on only during the equal exposure
        period, providing the same amount of exposure for all pixels in the ROI.

        Please see the camera specification for details on EEP mode.

        When enabled, EEP mode will be active or inactive depending on the exposure duration.
        Use :ref:`is_eep_supported<Ref-TLCamera-is_eep_supported>` to test whether the connected camera supports this
        mode.
        Use :ref:`eep_status<Ref-TLCamera-eep_status>` to see whether the mode is active, inactive, in bulb mode, or
        disabled.

        """
        try:
            eep_status = self.eep_status
            if eep_status == EEP_STATUS.DISABLED:
                return False
            else:
                return True
        except Exception as exception:
            _logger.error("Could not get is eep enabled; " + str(exception))
            raise exception

    @is_eep_enabled.setter
    def is_eep_enabled(self, is_eep_enabled):
        try:
            c_value = c_int(is_eep_enabled)
            error_code = self._sdk.tl_camera_set_is_eep_enabled(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_is_eep_enabled", error_code))
        except Exception as exception:
            _logger.error("Could not set is eep enabled; " + str(exception))
            raise exception

    @property
    def biny(self):
        """
        int:

        The current vertical binning value.

        """
        try:
            biny = c_int()
            error_code = self._sdk.tl_camera_get_biny(self._camera, biny)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_biny", error_code))
            return int(biny.value)
        except Exception as exception:
            _logger.error("Could not get bin y; " + str(exception))
            raise exception

    @biny.setter
    def biny(self, biny):
        try:
            c_value = c_int(biny)
            error_code = self._sdk.tl_camera_set_biny(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_biny", error_code))
        except Exception as exception:
            _logger.error("Could not set bin y; " + str(exception))
            raise exception

    @property
    def biny_range(self):
        """
        Range:

        The binning ratio in the Y direction can be determined with this property. By default, binning is set to one in
        both X and Y directions.
        This property is Read-Only. To set biny, see :ref:`biny<Ref-TLCamera-biny>`.

        """
        try:
            biny_min = c_int()
            biny_max = c_int()
            error_code = self._sdk.tl_camera_get_biny_range(self._camera, biny_min, biny_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_biny_range", error_code))
            return Range(int(biny_min.value), int(biny_max.value))
        except Exception as exception:
            _logger.error("Could not get biny range; " + str(exception))
            raise exception

    @property
    def gain(self):
        """
        int:

        Gain refers to the scaling of pixel values up or down for a given amount of light. This scaling is applied
        prior to digitization.
        To determine the valid range of values, use the :ref:`gain_range<Ref-TLCamera-gain_range>` property.
        If the :ref:`gain_range<Ref-TLCamera-gain_range>` maximum is zero, then Gain is not supported for the connected
        camera.
        The units of measure for Gain can vary by camera. Please consult the data sheet for the specific camera model.
        Query the :ref:`gain_range<Ref-TLCamera-gain_range>` property to determine the possible values.

        """
        try:
            gain = c_int()
            error_code = self._sdk.tl_camera_get_gain(self._camera, gain)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_gain", error_code))
            return int(gain.value)
        except Exception as exception:
            _logger.error("Could not get gain; " + str(exception))
            raise exception

    @gain.setter
    def gain(self, gain):
        try:
            c_value = c_int(gain)
            error_code = self._sdk.tl_camera_set_gain(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_gain", error_code))
        except Exception as exception:
            _logger.error("Could not set gain; " + str(exception))
            raise exception

    @property
    def black_level(self):
        """
        int:

        Black level adds an offset to pixel values. If the connected camera supports black level, the
        :ref:`black_level_range<Ref-TLCamera-black_level_range>` will have a maximum greater than zero.

        """
        try:
            black_level = c_int()
            error_code = self._sdk.tl_camera_get_black_level(self._camera, black_level)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_black_level", error_code))
            return int(black_level.value)
        except Exception as exception:
            _logger.error("Could not get black level; " + str(exception))
            raise exception

    @black_level.setter
    def black_level(self, black_level):
        try:
            c_value = c_int(black_level)
            error_code = self._sdk.tl_camera_set_black_level(self._camera, c_value)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_set_black_level", error_code))
        except Exception as exception:
            _logger.error("Could not set black level; " + str(exception))
            raise exception

    @property
    def black_level_range(self):
        """
        Range:

        Black level adds an offset to pixel values. If black level is supported by a camera model, then this property
        will have an upper range higher than zero.

        black_level_range indicates the available values that can be used for the
        :ref:`black_level<Ref-TLCamera-black_level>` property.

        """
        try:
            black_level_min = c_int()
            black_level_max = c_int()
            error_code = self._sdk.tl_camera_get_black_level_range(self._camera, black_level_min, black_level_max)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_black_level_range", error_code))
            return Range(int(black_level_min.value), int(black_level_max.value))
        except Exception as exception:
            _logger.error("Could not get black level range; " + str(exception))
            raise exception

    @property
    def image_width_pixels(self):
        """
        int:

        This property provides the image width in pixels. It is related to ROI width. This property is Read-Only.

        """
        try:
            image_width_pixels = c_int()
            error_code = self._sdk.tl_camera_get_image_width(self._camera, image_width_pixels)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_get_image_width", error_code))
            return int(image_width_pixels.value)
        except Exception as exception:
            _logger.error("Could not get image width; " + str(exception))
            raise exception

    @property
    def image_height_pixels(self):
        """
        int:

        This property provides the image height in pixels. It is related to ROI height. This property is Read-Only.

        """
        try:
            image_height_pixels = c_int()
            error_code = self._sdk.tl_camera_get_image_height(self._camera, image_height_pixels)
            if error_code != 0:
                raise TLCameraError(_create_c_failure_message(self._sdk, "tl_camera_image_height_pixels", error_code))
            return int(image_height_pixels.value)
        except Exception as exception:
            _logger.error("Could not get image height; " + str(exception))
            raise exception


""" Error handling """


class TLCameraError(Exception):
    def __init__(self, message):
        _logger.debug(message)
        super(TLCameraError, self).__init__(message)
