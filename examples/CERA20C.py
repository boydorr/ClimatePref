#!/usr/bin/env python
from ecmwfapi import ECMWFDataServer
server = ECMWFDataServer()

def retrieve_cera20c_edmm():
    """
       A function to demonstrate how to iterate efficiently over all months,
       for a list of years, for a CERA-20C synoptic monthly means request.
       You can extend the number of years to adapt the iteration to your needs.
       You can use the variable 'target' to organise the requested data in files as you wish."
    """
    yearStart = 1901
    yearEnd = 2010
    monthStart = 1
    monthEnd = 12

    months = range(monthStart, monthEnd + 1)
    years = range(yearStart, yearEnd+1)
    decades = list(set([divmod(i, 10)[0] for i in years]))
    decades = [x * 10 for x in decades]
    decades.sort()
    print 'Decades:', decades
    # loop through decades and create a month list
    for d in decades:
        requestDates=''
        for y in years:
            if ((divmod(y,10)[0])*10) == d:
                for m in months:
                    requestDates = requestDates+str(y)+(str(m)).zfill(2)+'01/'
        requestDates = requestDates[:-1]
        print 'Requesting dates: ', requestDates
        target = 'cera_20C_edmo_temp_%d'% (d)    # specifies the output file name
        print 'Output file: ', target
        era20c_edmm_sfc_request(requestDates, target)

    #requestMonthList = []
    #for year in list(range(yearStart, yearEnd + 1)):
    #    for month in list(range(monthStart, monthEnd + 1)):
    #        requestMonthList.append('%04d-%02d-01' % (year, month))
    #requestMonths = "/".join(requestMonthList)
    # we submit a data request for the current years
    #target = "cera20c_edmm_an_%s_sfc.grb" % (yearStart)
    #era20c_edmm_sfc_request(requestMonths, target)


    # edmm is arranged by years, so we iterate over the years
    #for year in list(range(yearStart, yearEnd + 1)):
    #    requestMonthList = []
    #    for month in list(range(monthStart, monthEnd + 1)):
    #        requestMonthList.append('%04d-%02d-01' % (year, month))
    #    requestMonths = "/".join(requestMonthList)
    #    # we submit a data request for the current year
    #    target = "cera20c_edmm_an_%s_sfc.grb" % (year)
    #    era20c_edmm_sfc_request(requestMonths, target)

def era20c_edmm_sfc_request(requestMonths, target):
    """
        A CERA-20C request for analysis, sfc data.
        You can change the keywords below to adapt it to your needs.
        (eg add or remove  parameters, times etc)
    """
    server.retrieve({
        "class": "ep",
        "dataset": "cera20c",
        "stream": "edmo",
        "expver": "1",
        "type": "an",
        "number": "0",
        "levtype": "sfc",
        "format": "netcdf",
        "grid": "0.375/0.375",
        "date": requestMonths,
        "param": "167.128",
        "target": target
    })

if __name__ == '__main__':
    retrieve_cera20c_edmm()
