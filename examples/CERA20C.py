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
    yearEnd = 1909
    monthStart = 1
    monthEnd = 12

    # edmm is arranged by years, so we iterate over the years
    for year in list(range(yearStart, yearEnd + 1)):
        requestMonthList = []
        for month in list(range(monthStart, monthEnd + 1)):
            requestMonthList.append('%04d-%02d-01' % (year, month))
        requestMonths = "/".join(requestMonthList)
        # we submit a data request for the current year
        target = "cera20c_edmm_an_%s_sfc.grb" % (year)
        era20c_edmm_sfc_request(requestMonths, target)

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
        "levtype": "sfc",
        "date": requestMonths,
        "param": "167.128",
        "target": target
    })

if __name__ == '__main__':
    retrieve_cera20c_edmm()
