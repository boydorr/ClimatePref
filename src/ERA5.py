import cdsapi

c = cdsapi.Client()
c.retrieve('reanalysis-era5-complete', {
    'class'   : 'ea',
    'expver'  : '1',
    'stream'  : 'moda',
    'type'    : 'an',
    'param'   : '167.128',
    'levtype' : 'sfc',
    'date'    : '2018-01-01',
    'decade'  : '2010',
}, 'monthly-mean-daily-mean-temp-an-sfc.grib')

def retrieve_era5():

    yearStart = 1979
    yearEnd = 2018
    months = [1,2,3,4,5,6,7,8,9,10,11,12]

    years = range(yearStart, yearEnd+1)
    print('Years: ',years)
    decades = list(set([divmod(i, 10)[0] for i in years]))
    decades = [x * 10 for x in decades]
    decades.sort()
    print('Decades:', decades)
    # loop through decades and create a month list
    for d in decades:
        requestDates=''
        for y in years:
            if ((divmod(y,10)[0])*10) == d:
                for m in months:
                    requestDates = requestDates+str(y)+(str(m)).zfill(2)+'01/'
        requestDates = requestDates[:-1]
        print('Requesting dates: ', requestDates)
        target = 'era_5_moda_%d'% (d)    # specifies the output file name
        print('Output file: ', target)
        era5_request(requestDates, d, target)

# the actual data request
def era5_request(requestDates, decade, target):
        c.retrieve('reanalysis-era5-complete', {
        "class": "ea",              # do not change
        'expver': '1',              # do not change
        'stream': 'moda',           # monthly means of daily means
        'type': 'an',               # analysis (versus forecast, fc)
        'levtype': 'sfc',           # surface data (versus pressure level, pl, and model level, ml)
        'param': '167.128',         # here: sea surface temperature (param 34) and mean sea level pressure (param 151)
        'format': 'netcdf',         # get output in netcdf; only works with regular grids; for GRIB remove this line
        'date': requestDates,       # dates, set automatically from above
        'decade': decade,           # decade set automatically from above
        'target': target            # output file name, set automatically from above
    })

if __name__ == '__main__':
    retrieve_era5()
