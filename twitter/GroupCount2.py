import numpy as np
from datetime import timedelta,datetime,timezone,tzinfo
from ReverseDatenum import zero_to_dates
import pytz 
import pandas as pd
from Datenum import days_from_zero


def serial_to_datetime(serial_date,timezone='Europe/London'):
    # Define the reference date as January 1, 0001
    ref_date = datetime(1, 1, 1)

 # Convert the serial date to a timedelta
    delta = timedelta(days=serial_date-365)  # Subtract one year to account for the year 0000

    # Add the timedelta to the reference date
    date = ref_date + delta

    # Make the datetime object timezone-aware
    date = date.replace(tzinfo=pytz.UTC)

    # Convert the date to the specified timezone
    date = date.astimezone(pytz.timezone(timezone))

    return date



def datetime_from_days(day_count: int, tz: tzinfo = timezone.utc) -> datetime:
    ref_date = datetime(1, 1, 1, tzinfo=tz)
    day_count_adjusted = day_count - 365
    return ref_date + timedelta(days=day_count_adjusted)





def groupcount(time, tweet, follower, timestep, del_t, start, stop, total_tweets):
    total_tweets = len(time)    # Correct length = 929 
    run_time = (stop - start) 
    run_time = run_time.total_seconds() / (24 * 60 * 60)  
    del_t = round(del_t, 4)
    del_t = timedelta(days=del_t)
    dataTime = [None]*(timestep+1)
    dataTime[0] = start-del_t
    dataTime[1] = start
    for i in range(2, timestep+1):
        dataTime[i] = dataTime[i-1] + del_t
    dataInf = np.zeros(timestep+1)
    dataSusc = np.zeros(timestep+1)
    timeupper = stop
    timelower = stop-del_t  
    time = time.dt.to_pydatetime()
    #time_zone = pytz.timezone('Europe/London')
    #time = [x.astimezone(time_zone) for x in time]

    for j in range(timestep+1, 1, -1):
        index = [(x >= timelower) & (x < timeupper) for x in time]
        for x in time[:20]:
            print(f"For x = {x}, (x >= timelower) & (x < timeupper) is {(x >= timelower) & (x < timeupper)}")
        temp_susc = 0
        temp_inf = 0
        for i in range(total_tweets): #time is not falling in the range between timelower and time upper
            if index[i]: #this for loop is not executing as it should, 
                temp_susc += follower[i]
                temp_inf += 1
        dataInf[j-1] = temp_inf
        dataSusc[j-1] = temp_susc
        timeupper = timelower
        timelower = timelower - del_t
# Next issue: dateInf and dateSusc are empty arrays, need to figure out why - fixed (manually)
    idx = np.nonzero(dataInf)[0]
    dataInf = dataInf[idx[0]-1:]
    dataSusc = dataSusc[idx[0]-1:]
    dataTime = dataTime[idx[0]-1:]
    dataInf[0] = 1
    dataSusc[0] = dataSusc[1]
    
    # Convert datetime back to timestamp
    dataTime = [item.to_pydatetime() for item in dataTime]

    return dataTime, dataSusc, dataInf