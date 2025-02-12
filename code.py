import time
import math
import board
import busio
#from adafruit_icm20x import ICM20948
from adafruit_lsm6ds.ism330dhcx import ISM330DHCX as ISM330
from adafruit_lis3mdl import LIS3MDL

# Initialize I2C bus and sensor
i2c = busio.I2C(board.SCL, board.SDA)
sensor_mag = LIS3MDL(i2c)
sensor_gyro = ISM330(i2c)

# Gyro default scale 250 dps. Convert to radians/sec subtract offsets
Gscale = (math.pi / 180.0) * 0.00763  # 250 dps scale sensitivity = 131 dps/LSB
G_offset = [-9.6516484e-05, -0.00660651, -0.00413769]#[74.3, 153.8, -5.5]


# Accel scale: divide by 16604.0 to normalize
A_B = [-0.07, -0.21, 0.04]
A_Ainv = [
    [103.77675 , -1.27451 , -1.51896],
    [-1.27451 , 101.14016 , 0.08487],
    [-1.51896 , 0.08487 , 102.24004]
]



# Mag scale divide by 369.4 to normalize
M_B = [-16.86035323, 17.56994745, -35.50285682]
M_Ainv = [
    [ 23.41701 , 0.80824 , -0.04391],
 [ 0.80824 , 24.26587 , 0.01451],
 [-0.04391 , 0.01451 , 23.46526]
]

# Local magnetic declination in degrees
declination = -0.13333#-14.84

# Mahony filter parameters
Kp = 50.0
Ki = 0.0

# Quaternion
q = [1.0, 0.0, 0.0, 0.0]
yaw, pitch, roll = 0.0, 0.0, 0.0

# Timers
last = time.monotonic()
PRINT_SPEED = 0.3  # seconds between angle prints
lastPrint = time.monotonic()

def vector_dot(a, b):
    return sum(x * y for x, y in zip(a, b))

def vector_normalize(a):
    mag = math.sqrt(vector_dot(a, a))
    return [x / mag for x in a]

def get_scaled_IMU():
    sensor_gyro.acceleration
    sensor_gyro.gyro
    sensor_mag.magnetic

    Gxyz = [
        Gscale * (sensor_gyro.gyro[0] - G_offset[0]),
        Gscale * (sensor_gyro.gyro[1] - G_offset[1]),
        Gscale * (sensor_gyro.gyro[2] - G_offset[2])
    ]

    Axyz = [sensor_gyro.acceleration[0], sensor_gyro.acceleration[1], sensor_gyro.acceleration[2]]
    Mxyz = [sensor_mag.magnetic[0], sensor_mag.magnetic[1], sensor_mag.magnetic[2]]

    # Apply accel offsets (bias) and scale factors from Magneto
    temp = [(Axyz[i] - A_B[i]) for i in range(3)]
    Axyz = [
        sum(A_Ainv[i][j] * temp[j] for j in range(3)) for i in range(3)
    ]
    Axyz = vector_normalize(Axyz)

    # Apply mag offsets (bias) and scale factors from Magneto
    temp = [(Mxyz[i] - M_B[i]) for i in range(3)]
    Mxyz = [
        sum(M_Ainv[i][j] * temp[j] for j in range(3)) for i in range(3)
    ]
    Mxyz = vector_normalize(Mxyz)

    return Gxyz, Axyz, Mxyz

def MahonyQuaternionUpdate(ax, ay, az, gx, gy, gz, mx, my, mz, deltat):
    global q
    eInt = [0.0, 0.0, 0.0]
    q1, q2, q3, q4 = q

    # Measured horizon vector = a x m (in body frame)
    hx = ay * mz - az * my
    hy = az * mx - ax * mz
    hz = ax * my - ay * mx
    norm = math.sqrt(hx * hx + hy * hy + hz * hz)
    if norm == 0.0:
        return
    norm = 1.0 / norm
    hx *= norm
    hy *= norm
    hz *= norm

    # Estimated direction of Up reference vector
    ux = 2.0 * (q2 * q4 - q1 * q3)
    uy = 2.0 * (q1 * q2 + q3 * q4)
    uz = q1 * q1 - q2 * q2 - q3 * q3 + q4 * q4

    # Estimated direction of horizon (West) reference vector
    wx = 2.0 * (q2 * q3 + q1 * q4)
    wy = q1 * q1 - q2 * q2 + q3 * q3 - q4 * q4
    wz = 2.0 * (q3 * q4 - q1 * q2)

    # Error is the summed cross products of estimated and measured directions of the reference vectors
    ex = (ay * uz - az * uy) + (hy * wz - hz * wy)
    ey = (az * ux - ax * uz) + (hz * wx - hx * wz)
    ez = (ax * uy - ay * ux) + (hx * wy - hy * wx)

    if Ki > 0.0:
        eInt[0] += ex
        eInt[1] += ey
        eInt[2] += ez
        gx += Ki * eInt[0]
        gy += Ki * eInt[1]
        gz += Ki * eInt[2]

    gx += Kp * ex
    gy += Kp * ey
    gz += Kp * ez

    gx *= 0.5 * deltat
    gy *= 0.5 * deltat
    gz *= 0.5 * deltat
    qa, qb, qc = q1, q2, q3
    q1 += (-qb * gx - qc * gy - q4 * gz)
    q2 += (qa * gx + qc * gz - q4 * gy)
    q3 += (qa * gy - qb * gz + q4 * gx)
    q4 += (qa * gz + qb * gy - qc * gx)

    norm = math.sqrt(q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4)
    norm = 1.0 / norm
    q = [q1 * norm, q2 * norm, q3 * norm, q4 * norm]

while True:
    Gxyz, Axyz, Mxyz = get_scaled_IMU()

    Mxyz[1] = -Mxyz[1]
    Mxyz[2] = -Mxyz[2]

    now = time.monotonic()
    deltat = now - last
    last = now

    MahonyQuaternionUpdate(Axyz[0], Axyz[1], Axyz[2], Gxyz[0], Gxyz[1], Gxyz[2], Mxyz[0], Mxyz[1], Mxyz[2], deltat)

    if time.monotonic() - lastPrint > PRINT_SPEED:
        roll = math.atan2((q[0] * q[1] + q[2] * q[3]), 0.5 - (q[1] * q[1] + q[2] * q[2]))
        pitch = math.asin(2.0 * (q[0] * q[2] - q[1] * q[3]))
        yaw = math.atan2((q[1] * q[2] + q[0] * q[3]), 0.5 - (q[2] * q[2] + q[3] * q[3]))

        yaw *= 180.0 / math.pi
        pitch *= 180.0 / math.pi
        roll *= 180.0 / math.pi
        roll += 45
        #yaw = -(yaw + declination)
        #if yaw < 0:
        #    yaw += 360.0
        #if yaw >= 360.0:
        #    yaw -= 360.0

        print(f"yaw: {yaw:.0f}, pitch: {pitch:.0f}, roll: {roll:.0f}")
        lastPrint = time.monotonic()