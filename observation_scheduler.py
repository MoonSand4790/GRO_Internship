import numpy as np
from datetime import datetime, timezone, timedelta, time
import pandas as pd
from zoneinfo import ZoneInfo

## Co-ordinate Changes
# Changing the Co-ordinate system from Galactic to Equatorial system
def gal_to_eq(l,b):
    # Required North Galactic Pole values
    l_ngp = np.radians(123.0)
    del_ngp = np.radians(27.133)
    ra_ngp = np.radians(192.85)
    
    l = np.radians(l)
    b = np.radians(b)

    # Declination
    decl_sin = np.cos(l - l_ngp) * np.cos(b) * np.cos(del_ngp) + np.sin(b) * np.sin(del_ngp)
    decl = np.arcsin(decl_sin)
    
    # Right Accession
    ra_sin = (np.sin(l_ngp - l)*np.cos(b))/np.cos(decl)
    ra_cos = (np.sin(decl)*np.cos(del_ngp) - np.cos(l_ngp - l)*np.cos(b))/(np.cos(decl)*np.sin(del_ngp))
    ra_s = np.arcsin(ra_sin) #+ ra_ngp
    
    # To determine the orientation in the sky 
    if ra_cos < 0:
        offset = 180 - np.degrees(ra_s)
    else:
        offset = np.degrees(ra_s)

    ra = (offset + np.degrees(ra_ngp)) % 360
    
    return ra, np.degrees(decl)

# Changing the Co-ordinate system from Equatorial to Horizontal system
def eq_to_hor(ra, dec, lat, longi, dt):
    # Local Sidereal Time
    lst = lst_calc(longi, dt)
    
    # Hour Angle (H is positive for West; Negative, East)
    h_angle = np.radians((lst - ra) % 360)
    
    dec_rad = np.radians(dec)
    lat_rad = np.radians(lat)

    # Altitude
    sin_alt = (np.sin(dec_rad) * np.sin(lat_rad) + 
               np.cos(dec_rad) * np.cos(lat_rad) * np.cos(h_angle))
    alt = np.arcsin(sin_alt)

    # Azimuth (South towards West)
    y = np.sin(h_angle)
    x = (np.cos(h_angle) * np.sin(lat_rad) - np.tan(dec_rad) * np.cos(lat_rad))
    az_rad = np.arctan2(y, x)
    
    # Shifting origin of azimuth to North
    az_deg = (np.degrees(az_rad) + 180) % 360

    return np.degrees(alt), az_deg

## Local Sidereal Time Calculation
# Calculation of Julian Day
def julian_day(year, month, day, hour, minute, sec):
    a = (14 - month) // 12
    y = year + 4800 - a
    m = month + 12 * a - 3
    
    jdn = day + (153 * m + 2) // 5 + 365 * y + y // 4 - y // 100 + y // 400 - 32045
    
    jd = jdn + (hour - 12) / 24 + minute / 1440 + sec / 86400
    return jd

# Calculation of Local Sidereal Time
def lst_calc(east_long, datetime):
    jd = julian_day(datetime.year, datetime.month, datetime.day, datetime.hour, datetime.minute, datetime.second)
    jdt = int(jd + 0.5)-0.5

    d = jdt - 2451545
    T = d / 36525
   
    gmst = 280.46061837 + 360.98564736629*(jd-2451545.0) + 0.000387933*(T**2) - (T**3)/38710000
    gmst = gmst % 360

    lmst = (gmst + east_long) % 360
    return lmst


## Main Code
start_time = datetime(2026,2,24,4,0,0,0)
start_time = start_time.astimezone(timezone.utc)
end_time = datetime(2026,2,24,9,0,0,0)
end_time = end_time.astimezone(timezone.utc)

# Latitude, Longitude and (assumed) Galactic latitude
b = 0.0
longi = 77.51
lat = 13.61

# Code runs for a set range of Galactic longitude (l) and saves for the Time range to a csv file
for l in range(18,90,6):
    rows = []
    current_time = start_time
    while current_time < end_time:
        ra, dec = gal_to_eq(l,b)
        alt, az = eq_to_hor(ra,dec,lat,longi,current_time)
        rows.append({
            "date":current_time.astimezone(ZoneInfo("Asia/Kolkata")).date(),
            "time":current_time.astimezone(ZoneInfo("Asia/Kolkata")).time(),
            "ra":ra,
            "dec":dec,
            "alt":alt,
            "az":az
        })
        current_time += timedelta(minutes=15)

    df = pd.DataFrame(rows)
    df.to_csv(f"obs_l{l}_date_month.csv", index=False)

print('Done')
