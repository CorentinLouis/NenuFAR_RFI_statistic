#−∗− coding : utf−8 −∗−

import datetime
import numpy
import pytz
from astropy.time import Time, TimeDelta

@numpy.vectorize
def doy_to_ymd(doy,hours,minutes,second):
    ### Function to change a (doy (in yyyyddd format), hours, minutes, second) to a datetime object datetime.datetime(YYYY, mm, dd, HH, MM, SS)
    return datetime.datetime.strptime(doy+hours+minutes+second,'%Y%j%H%M%S.%f')


def doy_to_datetime(time_doy):
    ### Function to change a doy (in floating yyyyddd format) to a datetime object datetime.datetime(YYYY, mm, dd, HH, MM, SS)
    time_hours = [int((itime-int(itime))*24) for itime in (time_doy)]
    time_minutes = [int(((time_doy[itime]-int(time_doy[itime]))*24-time_hours[itime])*60) for itime in range(len(time_doy))]
    time_seconds = [int((((time_doy[itime]-int(time_doy[itime]))*24-time_hours[itime])*60-time_minutes[itime])*60) for itime in range(len(time_doy))]
    time = [datetime.datetime.strptime(f'{int(time_doy[itime])}T{time_hours[itime]:02d}:{time_minutes[itime]:02d}:{time_seconds[itime]:02d}', "%Y%jT%H:%M:%S") for itime in range(len(time_doy))]
    return(time)


@numpy.vectorize
def doy_float_to_ymd(doy,hours):
    ### Function to transform a doy YYYYDDD in floating format to a datetime object
    ymd = datetime.datetime.strptime((doy)[0:7],'%Y%j')
    daydec = datetime.timedelta(hours=(hours))
    return (ymd + daydec)

@numpy.vectorize
def return_hour_minute(date):
    ### Function to return only the hour and minute of a datetime object
    result = datetime.datetime.strftime(date,'%H:%M')
    return (result)

@numpy.vectorize
def datetime_to_timestamp(datetime_table):
    ### Function to return time in floating format (from a datetime object)
    return Time(datetime_table, format="datetime").unix

@numpy.vectorize
def timestamp_to_datetime(timestamp_table):
    ### Function to return time in datetime format (from a timestamp object)
    result = Time(timestamp_table, format="unix").datetime
    return (result)

def doy_specific_year_to_yyyyddd(doy, origin):  
### Function to change "day of a specific year" format to yyyyddd ###

    aa = numpy.arange(61, dtype=float)+origin  # array of years starting from year of origin
    deb = numpy.zeros([61], dtype=float)  # zeros
    for i in range(1, len(deb)):  # categorising start point for each year
        if i % 4 == 1:
            deb[i:] = deb[i:]+366.
        else:
            deb[i:] = deb[i:]+365.

    yyyyddd = numpy.zeros(len(doy), dtype=float)

    for i in range(0, len(doy)):
        j = doy[i]-deb
        yyyyddd[i] = (aa[j >= 1][-1])*1000.+j[j >= 1][-1]

    return(yyyyddd)

@numpy.vectorize
def from_utc_to_local_time(time_datetime, timezone = "Europe/Paris"):
    """
    Transfrom the UTC time to a local (timezone) time
    Availabel timezone list accessible through: pytz.all_timezones
    """
    tz = pytz.timezone(timezone)
    local_time_datetime = tz.fromutc(time_datetime)
    return(local_time_datetime)


@numpy.vectorize
def datetime_to_float_since1970(time_datetime):
    """
    Return the time in days since 1970-01-01T00:00:00 (timestamp_0)
    """
    epoch = datetime.datetime.utcfromtimestamp(0)
    time_float_days =  (time_datetime.replace(tzinfo=None) - epoch).total_seconds()/60/60/24
    return time_float_days

