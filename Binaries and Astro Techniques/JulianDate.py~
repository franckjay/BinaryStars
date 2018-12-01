#########################################################################
##Jay Franck
##AST 680 Homework 3a
##
##Introduction: This quick python program is designed to easily convert
##Gregorian Dates in the format: YYYY MMM DD into the correct Julian Date.
##The Julian Date is especially useful in all aspects of computer-based
##systems as it is completely numerical and does not rely on leap year rules,
##months of the year, etc. It is useful especially in the realm of astronomy
##for binary/periodic systems(finding phases, predicting eclipses), as the 
##computer modelling for these systems are largely based on time. This program
##has been written so that it can easily be ported to other input styles or
##output requirements.
#########################################################################
import math

def Month(month):
    '''Input into this subroutine is the first three letters of the month'''
    ##Output will be the month in integer format i.e. March as 03
    ##By using the command 'in' the program is still able to differentiate
    ##the month even if input is not in three letter format(i.e. Mar for March).
    ##Obviously there are cleaner and more efficient ways of computing the month
    ##but this is the quick and dirty way without dealing with Python's
    ##dictionary system.

    Mo=0
    if 'Jan' in month:
        Mo=1
    elif 'Feb' in month:
        Mo=2
    elif 'Mar' in month:
        Mo=3
    elif 'Apr' in month:
        Mo=4
    elif 'May' in month:
        Mo=5
    elif 'Jun' in month:
        Mo=6
    elif 'Jul' in month:
        Mo=7
    elif 'Aug' in month:
        Mo=8
    elif 'Sep' in month:
        Mo=9
    elif 'Oct' in month:
        Mo=10
    elif 'Nov' in month:
        Mo=11
    elif 'Dec' in month:
        Mo=12
    else:
        print "Incorrect input. Please run program again"
    return Mo

def JulianDate(year,month,day):
    ##This subroutine will take the Gregorian Date in YYYY MonthName DD and
    ##convert it to the Julian Date using the subroutine from above. This will
    ##be using the algorithm by R.G. Tantzen(Comm. ACM, vol. 6, #8; adapted by
    ##George Diedrich in 1964 for use by amateur astronomers).
    year = int(year)
    day = int(day)
    Mo= Month(month)
    ##Significance of 'A': Since the Julian Date is a base ten system(decimal)
    ##the month terms needs to be converted from a 12 number system into a
    ##decimal system to be utilized efficiently in this algorithm.
    A=int(Mo-3)
    if A < 0:
        B = A+12
        C=int(year-1)
    else:
        B=A
        C=year
    D = int(C/100)
    E = int(C-(100*D))
    F = int((146097*D)/4)
    G = int((1461*E)/4)
    H = int(((153*B)+2)/5)
    ##Finally, for the Julian Date Calculation!!
    JDoo = F + G + H +1721119 + day - 0.5
    ##No correction for JDut decimal days, as each of the days in this program
    ##are at 00 hours UT. However, if needed, the implementation of this would
    ##be a trivial matter. The program would be modified to perform (days-days%1)
    ##which I believe is the mod function in other languages. This would get the
    ##remainder days IFF the input was in UT time. The program would then
    ##multiply the remainder days by 24 in order to get the decimal hours
    ##and then would implement the formula JDut=JDoo+(UT/24)
    output= "%4i %4s %5i %3i %4i %6i %3i %3i %i %4i %5i %.1f"%(year,month,day,A,B,C,D,E,F,G,H,JDoo)
    print output
    return JDoo

##This is the main body of the program. It features a single array with each
##of the dates as a string in the format provided. It would be very easy to
##alter this input to suit a variety of program requirements, including manual
##input and input from a file.
GD=['1999 Dec 31','2000 Jan 1','2000 Mar 1','2010 Sep 30','2100 Mar 1','1582 Oct 4','1582 Oct 15']

##This for loop iterates through each of the values(i), splits the individual
##pieces of data into a new temporary array (YYYYMMMDD) and runs it through
##the JulianDate subroutine. It then prints the Gregorian Date and Julian Date.
header="%4s %4s %4s %4s %4s %4s %5s %3s %4s %4s %8s %6s"%("Year","Mo","Day","A","B","C","D","E","F","G","H","JDoo")
print header
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
for i in GD:
    YYYYMMMDD=i.split()
    JulianDate(YYYYMMMDD[0],YYYYMMMDD[1],YYYYMMMDD[2])



###############################################################################
###############################################################################
## Questions at the end of the assignment:
##
## a)Comparing my answers with those from the US Naval Observatory, I found
## that my answers were all correct EXCEPT for on October 4, 1582, where my
## answer was exactly then days short of the USNO date. This can be accounted
## for as the USNO calculator takes into account the 10 "lost" days between
## October 4 and October 15, while my program does not (although it could be
## adjusted to easily).
##
## b) The Julian Calendar was the calendar adopted by Julius Caesar which had
## an average year of 365.25 days. Thus, every four years would be a leap year
## with an extra day (Fraknoi et al.71). However, this approximation led to
## a difference in the day of spring by about 10 days in the time of Pope
## Gregory. The Julian Date, on the other hand, is simply a decimal based date
## based on the number of days since the assumed start of 'recorded' history.
##
## c)The Gregorian calendar was introduced in order to provide a more accurate
## calendar than the former Julian Calendar, as it was based on more precise
## values of the length of the year, closer to the tropical year (365.2425 days)
## instead of the Julian Calendar year of 365.25 days (Fraknoi et al. 71).
##
## References: "Voyages to the Stars and Galaxies." Fraknoi,Morrison,Wolff.
##  Copyright 2006, Thomson,Brook,Cole Publishing.
##
## USNO Calendar Date to Julian Date Converter. http://www.usno.navy.mil/USNO/astronomical-
##  applications/data-services/cal-to-jd-conv. Accessed 25 September 2010.
