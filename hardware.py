import serial
import time
import os

os.environ['PATH'] = 'C:\\Program Files\\Thorlabs\\APT\\APT Server;' + os.environ['PATH']

# print(f'> loading hardware.py ...')

usb_oriel = 'COM10'
oriel_motor = 'M1'
usb_arduino = 'COM8'
usb_ds345 = 'COM7'
usb_mdt694b = 'COM9'
usb_kbd101 = ''

t_wait = 0.5


class Arduino:
    def __init__(self):
        self.usb = usb_arduino
        self.port = serial.Serial(self.usb, baudrate=9600, timeout=1)
        self.led = 0
        self.tcam = 0
        print(f'> hw.Arduino.init: usb @ {self.usb}')

    def write(self, cmd, val):
        cmd = cmd.upper()
        msg = cmd.encode() + bytes([val//256, val%256])
        self.port.write(msg)
        # print(f'> hw.Arduino.write: msg = {msg.decode().rstrip()}')

    def set_LT(self, LED=1, Tcam=100):
        time.sleep(t_wait)
        self.write('L', LED)
        time.sleep(t_wait)
        self.write('T', Tcam)
        time.sleep(t_wait)
        self.led = LED
        self.tcam = Tcam
        print(f'> hw.Arduino.set_LT: LED = {self.led}; Tcam = {self.tcam}')

    def set_LED(self, LED=0):
        '''LED_and_CamTrig.ino is set up as LED=1: pin 9; LED=2: pin 10; LED=3: pin 11'''
        time.sleep(t_wait)
        self.write('L', LED)
        time.sleep(t_wait)
        print(f'> hw.Arduino.set_LED: LED = {LED}')

    def pls_LED(self, LED=0, t_pls=100):
        '''t_pls in ms'''
        time.sleep(t_wait)
        self.write('L', LED)
        time.sleep(t_pls*1e-3)
        self.write('L', 0)
        time.sleep(t_pls*1e-3)
        print(f'> hw.Arduino.pls_LED: LED = {LED}')

    def set_Tcam(self, Tcam=100):
        time.sleep(t_wait)
        self.write('T', Tcam)
        time.sleep(t_wait)
        self.tcam = Tcam
        print(f'> hw.Arduino.set_Tcam: Tcam = {self.tcam}')

    def test_arduino(self):
        self.set_LT()
        self.set_Tcam()
        for i in range(3):
            self.pls_LED(LED=1, t_pls=500)


class MDT694B:
    def __init__(self):
        self.usb = usb_mdt694b
        self.port = serial.Serial(self.usb, baudrate=115200, timeout=1)
        self.port.flushInput()
        print(f'> hw.MDT694B.init: usb @ {self.usb}')

    def set_vpz(self, v_piezo, prn=True):
        msg = f'xvoltage={v_piezo:.3f}\r\n'.encode()
        self.port.write(msg)
        if prn:
            print(f'> hw.MDT694B.set_vpz: v_piezo = {v_piezo:.3f}; msg = {msg.decode().rstrip()}')
        time.sleep(t_wait)

    def get_vpz(self):
        msg = f'xvoltage?\r\n'.encode()
        self.port.write(msg)
        # print(f'msg = {msg}')
        v_piezo = self.port.readline()
        # print(f'v_piezo_1 = {v_piezo}')
        v_piezo = v_piezo[2:-3].decode()
        # print(f'v_piezo_2 = {v_piezo}')
        v_piezo = float(v_piezo)
        print(f'> hw.MDT694B.get_vpz: v_piezo = {v_piezo};  msg = {msg.decode().rstrip()}')
        return v_piezo

    def test_mdt(self):
        v_piezo = self.get_vpz()
        print(f'> test_mdt: V_piezo = {v_piezo}')
        v_req = input('> Enter V_req: ')
        self.set_vpz(v_req)
        v_piezo = self.get_vpz()
        print(f'> V_req = {v_req}; V_piezo = {v_piezo}')


class Quantalux:
    def __init__(self):
        from qcam_setup.thorlabs_tsi_sdk.tl_camera import TLCameraSDK
        from qcam_setup.windows_setup import configure_path
        import numpy as np

        configure_path()
        self.sdk = TLCameraSDK()
        self.camera = self.sdk.open_camera(self.sdk.discover_available_cameras()[0])
        self.nxy = (1920, 1080)
        self.nx, self.ny = self.nxy
        self.cc = np.empty((self.ny, self.nx), dtype=np.uint16)
        self.bit = 16
        self.iframe = 0
        self.tframe = time.time()
        self.roi = (0, 0, 0, 0)
        self.rnxy = (0, 0)
        self.rnx, self.rny = self.rnxy
        print(f'> hw.Quantalux.init: nxy = {self.nxy}')

    def setup(self, t_exp=1.0, trig=1, roi=(0, 0, 0, 0)):
        from qcam_setup.thorlabs_tsi_sdk.tl_camera_enums import OPERATION_MODE, TRIGGER_POLARITY

        self.camera.disarm()
        self.camera.exposure_time_us = int(t_exp*1000)
        if trig:
            self.camera.frames_per_trigger_zero_for_unlimited = 1
            self.camera.operation_mode = OPERATION_MODE.HARDWARE_TRIGGERED
            self.camera.trigger_polarity = TRIGGER_POLARITY.ACTIVE_HIGH
            # print('> hw.Quantalux.setup: trig is true')
        else:
            self.camera.frames_per_trigger_zero_for_unlimited = 0
            # self.camera.operation_mode = OPERATION_MODE.SOFTWARE_TRIGGERED
            # print('> hw.Quantalux.setup: trig is none')

        if roi != (0, 0, 0, 0):
            # roi = (0, 0, self.nx, self.ny)
            # else:
            self.camera.roi = roi
            x0, y0, x1, y1 = (self.camera.roi.upper_left_x_pixels, self.camera.roi.upper_left_y_pixels,
                              self.camera.roi.lower_right_x_pixels, self.camera.roi.lower_right_y_pixels)
            self.rnxy = (x1-x0+1, y1-y0+1)
            self.rnx, self.rny = self.rnxy
            self.roi = (x0, y0, self.rnx, self.rny)

        self.camera.arm(frames_to_buffer=1)
        print(f'> hw.Quantalux.setup: t_exp = {t_exp}, trig = {trig}, roi = {self.roi}')

    def acquire_cc(self):
        nfrm = 0
        while True:
            frame = self.camera.get_pending_frame_or_null()
            nfrm += 1
            # print(f'> hw.Quantalux.acquire_cc: looping: frame = {frame}')
            if frame:
                # print(f'> hw.Quantalux.acquire_cc: got frame @ nfrm = {nfrm}: frame = {frame}')
                break
            if nfrm > 100:
                print(f'> hw.Quantalux.acquire_cc: nfrm = {nfrm} and still no frame ...')
                break
        self.cc = frame.image_buffer
        self.iframe = frame.frame_count
        self.tframe = time.time()
        # print(f'> hw.Quantalux.acquire_cc: iframe = {self.iframe}, tframe = {self.tframe}')
        # gf.what_is('cc', self.cc)

    def dispose(self):
        self.camera.disarm()
        self.camera.dispose()
        self.sdk.dispose()
        print(f'> hw.Quantalux.dispose: disposed ...')

    def test_qcam(self):
        import pycubelib.plotting_functions as pf

        self.setup(t_exp=1.0)
        self.acquire_cc()
        pf.plotAA(self.cc, sxy=(.5, .5))
        print(f'< cc = {self.cc[0, 0:20]}')
        self.dispose()


class DDR25:
    def __init__(self):
        import thorlabs_apt as apt
        # _, self.sn = apt.list_available_devices()[0]
        self.sn = 28250272
        self.motor = apt.Motor(self.sn)
        print(f'> hw.DDR25.init: SN = {self.sn}')

    def position(self):
        q = self.motor.position
        print(f'> DDR25: Q = {q:.4f}')
        return q

    def direction(self, d='+'):
        if d == '+':
            self.motor.move_velocity(1)
        else:
            self.motor.move_velocity(2)

    def goto(self, q):
        q0 = self.position()
        if q > q0:
            self.direction('+')
            print(f'> DDR25: goto {q:.4f} +')
        else:
            self.direction('-')
            print(f'> DDR25: goto {q:.4f} -')
        self.motor.move_to(q)
        time.sleep(.5)
        self.position()

    def test_ddr(self):
        print(f'< DDR25 SN = {self.sn}')
        q = self.position()
        print(f'< current position = {q}')
        self.direction()
        self.goto(q + 1.23)
        print(f'< current position = {self.position()}')


class Oriel:
    def __init__(self):
        self.usb = usb_oriel
        self.motor = oriel_motor
        self.port = serial.Serial(self.usb, timeout=1)
        self.velocity = 0.0
        self.zpos = 0.0
        if not self.port.is_open:
            self.port.open()
        self.port.write(b'R\r\n')
        msg = (self.motor+'\r\n').encode()
        self.port.write(msg)
        time.sleep(t_wait)
        print(f'> hw.Oriel.init: usb @ {self.usb}')

    def position(self):
        self.port.reset_input_buffer()
        self.port.write(b'A\r\n')
        z = self.port.read_until(size=10).decode()
        z = z[:z.find('\\')-1]
        try:
            self.zpos = float(z)
        except ValueError:
            self.zpos = 0.0
        print(f'> hw.Oriel.position: zpos = {self.zpos}')
        return self.zpos

    def speed(self, vel):
        time.sleep(t_wait)
        self.port.write(f'V{vel}\r\n'.encode())
        self.velocity = vel
        time.sleep(t_wait)
        print(f'> hw.Oriel.speed: velocity = {self.velocity}')

    def goto(self, z):
        time.sleep(t_wait)
        self.port.write(f'G{z}\r\n'.encode())
        time.sleep(t_wait)
        print(f'> hw.Oriel.goto: z = {z}')

    def stop(self):
        time.sleep(t_wait)
        self.port.write(b'S\r\n')
        time.sleep(t_wait)
        print(f'> he.Oriel.stop: stopped ...')

    def close(self):
        time.sleep(t_wait)
        self.stop()
        self.port.close()
        print(f'> hw.Oriel.close: closed ...')


class DS345:
    def __init__(self):
        self.usb = usb_ds345
        self.port = serial.Serial(self.usb, timeout=1)
        msg = 'AMPL 0\r'.encode()
        self.port.write(msg)
        # msg = 'OFFS 0.0\r'.encode()
        # self.port.write(msg)
        # msg = 'FREQ 0.1\r'.encode()
        # aaa = self.port.write(msg)
        aaa = self.port.write('*IDN?\r'.encode())
        bbb = self.port.readline()#.decode()
        print(f'> hw.DS345.init: usb @ {self.usb} ...')


if __name__ == '__main__':

    # Arduino().test_arduino()
    # MDT694B().test_mdt()
    # Quantalux().test_qcam()
    DDR25().test_ddr()

    pass






