from ecmwfapi import ECMWFDataServer

server = ECMWFDataServer()


def retrieve_era_interim(param, from_year, to_year, filename='era_interim', **kwargs):
    """
    A function to retrieve ERA-Interim monthly means data.
    For data retrieval efficiency, all data access is made per decade.
    Requests can be made for a particular parameter, `param`, across a
    year range, `from_year`:`to_year`, with additional request arguments,
    including:

    :param filename: the filename (will be suffixed with decade)
    :param stream : monthly means of daily means (moda) or
        monthly means of daily forecast accumulations (mdfa)
    :param type : analyis (an) versus forecast (fc)
    :param levtype : susurface data (sfc) versus pressure level (pl)
        and model level (ml)
    :param grid : horizontal resolution of output in degrees lat/lon
    :param format : netcdf file or GRIB
    :param step : time step of forecast model (fc only)
    """
    months = range(1, 13)
    years = range(from_year, to_year + 1)
    decades = sorted({round(y, -1) for y in years})

    # Loop through decades and create a request list for all months/years
    for d in decades:
        # Filter for years within the decade
        years_in_decade = list(filter(lambda y: round(y, -1) == d, years))

        # Set up all request months per decade in correct format
        request_dates = "/".join([f'{y}{m:02}01' for y in years_in_decade for m in months])

        # Create target file
        target = f'{filename}_{d}'
        print(f'Years: {years_in_decade}\nOutput file: {target}')
        era_interim_request(param, request_dates, d, target, **kwargs)


def era_interim_request(
        param, request_dates, decade, target,
        stream='moda', type='an', levtype='sfc',
        grid='0.75/0.75', format='netcdf', step='0-12'
):
    """
    Request for ERA-interim data through MARS.
    """
    server.retrieve({
        "class":   "ei",
        "dataset": "interim",
        'expver':  '1',
        'stream':  stream,
        'type':    type,
        'levtype': levtype,
        'param':   param,
        'grid':    grid,
        'format':  format,
        'date':    request_dates,
        'decade':  decade,
        'target':  target,
        'step':    step
    })
