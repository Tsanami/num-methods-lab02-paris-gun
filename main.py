import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
import math

# ====================== КОНСТАНТЫ И МОДЕЛИ ======================
class EarthModel:
    # Параметры Земли (WGS84)
    a = 6378136.3  # Экваториальный радиус [м]
    f = 1 / 298.257223563
    b = a * (1 - f)  # Полярный радиус [м]
    eq_radius = 6378136.3  # Экваториальный радиус [м]
    polar_radius = 6356752.3  # Полярный радиус [м]
    J2 = 1.0826359e-3  # Коэффициент несферичности
    GM = 3.986004418e14  # Гравитационный параметр [м^3/с^2]
    omega = 7.292115e-5  # Угловая скорость вращения Земли [рад/с]

class AtmosphereModel:
    @staticmethod
    def density(h):
        """Плотность воздуха в зависимости от высоты [кг/м^3]"""
        if h < 11000:
            return 1.225 * np.exp(-h / 8500)
        elif h < 25000:
            return 0.364 * np.exp(-(h - 11000) / 6700)
        else:
            return 0.088 * np.exp(-(h - 25000) / 9000)

    @staticmethod
    def temperature(h):
        """Температура воздуха в зависимости от высоты [K]"""
        if h < 11000:
            return 288.15 - 0.0065 * h
        elif h < 25000:
            return 216.65
        else:
            return 216.65 + 0.003 * (h - 25000)

    @staticmethod
    def speed_of_sound(h):
        """Скорость звука в зависимости от высоты [м/с]"""
        gamma = 1.4
        R = 287.0  # Газовая постоянная для воздуха [Дж/(кг·K)]
        return math.sqrt(gamma * R * AtmosphereModel.temperature(h))

class Projectile:
    def __init__(self, caliber, mass):
        self.caliber = caliber  # Калибр снаряда [м]
        self.mass = mass  # Масса снаряда [кг]
        self.area = math.pi * (caliber / 2) ** 2  # Площадь поперечного сечения [м^2]
        self.drag_coeff = 0.25  # Коэффициент сопротивления (упрощенная модель)

# ====================== ПРЕОБРАЗОВАНИЕ КООРДИНАТ ======================
def geodetic_to_geocentric(lon, lat, h):
    """
    Преобразование геодезических координат в геоцентрические
    :param lon: Долгота [рад]
    :param lat: Широта [рад]
    :param h: Высота [м]
    :return: (x, y, z) в метрах
    """
    a, b = EarthModel.a, EarthModel.b
    e_sq = 1 - (b**2 / a**2)
    N = a / math.sqrt(1 - e_sq * math.sin(lat)**2)
    
    x = (N + h) * math.cos(lat) * math.cos(lon)
    y = (N + h) * math.cos(lat) * math.sin(lon)
    z = (N * (1 - e_sq) + h) * math.sin(lat)
    
    return x, y, z

def geocentric_to_geodetic(x, y, z):
    """
    Преобразование геоцентрических координат в геодезические (итеративный метод)
    :param x, y, z: Координаты [м]
    :return: (lon, lat, h) в радианах и метрах
    """
    a, b = EarthModel.a, EarthModel.b
    e_sq = 1 - (b**2 / a**2)
    lon = math.atan2(y, x)
    
    # Начальное приближение
    p = math.sqrt(x**2 + y**2)
    lat = math.atan2(z, p * (1 - e_sq))
    
    # Итерационный процесс
    max_iter = 10
    tolerance = 1e-6
    for _ in range(max_iter):
        N = a / math.sqrt(1 - e_sq * math.sin(lat)**2)
        h = p / math.cos(lat) - N
        lat_new = math.atan2(z, p * (1 - e_sq * N / (N + h)))
        
        if abs(lat_new - lat) < tolerance:
            break
        lat = lat_new
    
    return lon, lat, h

# ====================== ФИЗИЧЕСКИЕ МОДЕЛИ ======================
def gravity_acceleration(x, y, z):
    """
    Гравитационное ускорение с учетом J2
    :return: (gx, gy, gz) [м/с^2]
    """
    r = math.sqrt(x**2 + y**2 + z**2)
    r_sq = r**2
    z_sq = z**2
    
    # Сферическая составляющая
    common = -EarthModel.GM / r**3
    
    # Поправка J2
    J2_term = 1.5 * EarthModel.J2 * (EarthModel.eq_radius**2 / r_sq)
    z_factor = (5 * z_sq / r_sq - 1)
    
    gx = common * x * (1 - J2_term * z_factor)
    gy = common * y * (1 - J2_term * z_factor)
    gz = common * z * (1 - J2_term * (5 * z_sq / r_sq - 3))
    
    return gx, gy, gz

def drag_acceleration(state, projectile, atmosphere):
    """
    Ускорение от аэродинамического сопротивления
    :param state: [x, y, z, vx, vy, vz]
    :return: (ax, ay, az) [м/с^2]
    """
    x, y, z, vx, vy, vz = state
    _, _, h = geocentric_to_geodetic(x, y, z)
    
    # Атмосферные параметры
    rho = atmosphere.density(h)
    a_sound = atmosphere.speed_of_sound(h)
    
    # Аэродинамические параметры
    v = math.sqrt(vx**2 + vy**2 + vz**2)
    mach = v / a_sound
    drag_force = 0.5 * rho * v**2 * projectile.area * projectile.drag_coeff
    
    # Направление противоположно скорости
    ax = -drag_force * vx / (projectile.mass * v) if v > 0 else 0
    ay = -drag_force * vy / (projectile.mass * v) if v > 0 else 0
    az = -drag_force * vz / (projectile.mass * v) if v > 0 else 0
    
    return ax, ay, az

# ====================== ДИНАМИКА СНАРЯДА ======================
def projectile_dynamics(t, state, projectile, atmosphere):
    """
    Правые части дифференциальных уравнений движения снаряда
    :param state: [x, y, z, vx, vy, vz]
    :return: Производные состояния
    """
    x, y, z, vx, vy, vz = state
    
    # Гравитационное ускорение
    gx, gy, gz = gravity_acceleration(x, y, z)
    
    # Аэродинамическое ускорение
    ax_drag, ay_drag, az_drag = drag_acceleration(state, projectile, atmosphere)
    
    # Суммарное ускорение
    ax_total = gx + ax_drag
    ay_total = gy + ay_drag
    az_total = gz + az_drag
    
    return [vx, vy, vz, ax_total, ay_total, az_total]

def calculate_range(elevation, init_conditions, projectile):
    """
    Расчет дальности полета для заданного угла возвышения
    :param elevation: Угол возвышения [градусы]
    :return: Дальность [м]
    """
    # Распаковка начальных условий
    lon0, lat0, h0, v0, azimuth = init_conditions
    elevation_rad = math.radians(elevation)
    azimuth_rad = math.radians(azimuth)
    
    # Начальная скорость в локальной системе
    vx_local = v0 * math.cos(elevation_rad) * math.sin(azimuth_rad)
    vy_local = v0 * math.cos(elevation_rad) * math.cos(azimuth_rad)
    vz_local = v0 * math.sin(elevation_rad)
    
    # Преобразование в геоцентрические координаты
    x0, y0, z0 = geodetic_to_geocentric(lon0, lat0, h0)
    
    # Матрица преобразования
    clon, slon = math.cos(lon0), math.sin(lon0)
    clat, slat = math.cos(lat0), math.sin(lat0)
    
    rotation_matrix = np.array([
        [-slon, -slat*clon, clat*clon],
        [ clon, -slat*slon, clat*slon],
        [    0,       clat,      slat]
    ])
    
    # Начальная скорость в глобальной системе
    vx0, vy0, vz0 = rotation_matrix @ np.array([vx_local, vy_local, vz_local])
    
    # Начальное состояние
    initial_state = [x0, y0, z0, vx0, vy0, vz0]
    atmosphere = AtmosphereModel()
    
    # Интегрирование траектории
    def event(t, state):
        _, _, h = geocentric_to_geodetic(state[0], state[1], state[2])
        return h  # Событие при h=0
    
    event.terminal = True
    event.direction = -1
    
    sol = solve_ivp(
        fun=lambda t, y: projectile_dynamics(t, y, projectile, atmosphere),
        t_span=[0, 1000],
        y0=initial_state,
        events=event,
        method='RK45',
        rtol=1e-6,
        atol=1e-9
    )
    
    # Расчет дальности
    if sol.t_events[0].size > 0:
        x_impact, y_impact, z_impact = sol.y[0:3, -1]
        lon_impact, lat_impact, _ = geocentric_to_geodetic(x_impact, y_impact, z_impact)
        
        # Вычисление расстояния по дуге большого круга
        dlat = lat_impact - lat0
        dlon = lon_impact - lon0
        a = math.sin(dlat/2)**2 + math.cos(lat0) * math.cos(lat_impact) * math.sin(dlon/2)**2
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        return EarthModel.a * c
    return 0

# ====================== ОПТИМИЗАЦИЯ И ЗАПУСК ======================
def optimize_elevation(init_conditions, projectile):
    """
    Поиск оптимального угла возвышения для максимальной дальности
    :return: Оптимальный угол [градусы], Максимальная дальность [м]
    """
    def objective(elevation):
        return -calculate_range(elevation, init_conditions, projectile)
    
    result = minimize_scalar(
        objective,
        bounds=(30, 70),
        method='bounded',
        options={'xatol': 0.1}
    )
    
    return result.x, -result.fun

# ====================== ПАРАМЕТРЫ МОДЕЛИРОВАНИЯ ======================
if __name__ == "__main__":
    # Параметры снаряда (Парижская пушка)
    projectile = Projectile(caliber=0.21, mass=94)  # 210 мм, 94 кг
    
    # Начальные условия (пример)
    init_conditions = (
        math.radians(8.0),   # Долгота старта [рад] (Германия)
        math.radians(50.0),  # Широта старта [рад]
        0.0,                 # Высота старта [м]
        1600,                # Начальная скорость [м/с]
        225                  # Азимут (Юго-Запад) [град]
    )
    
    # Оптимизация угла возвышения
    opt_elevation, max_range = optimize_elevation(init_conditions, projectile)
    
    print(f"Оптимальный угол возвышения: {opt_elevation:.2f}°")
    print(f"Максимальная дальность: {max_range/1000:.2f} км")
    
    # Дополнительная проверка для 55°
    range_55 = calculate_range(55, init_conditions, projectile)
    print(f"Дальность при 55°: {range_55/1000:.2f} км")