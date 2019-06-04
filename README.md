# mjdemetra
Matlab function to perform seasonal adjustment with JDemetra+

## Installation

- Make sure Matlab uses the appropiate Java version
(type ```version -java``` in Matlab to find out which version is used)
Matlab should use a Java SE version that is [compatible with
JDemetra+](https://github.com/jdemetra/jdemetra-app) 

- Edit the ```classpath.txt``` file (type ```which classpath.txt``` in Matlab to find its location)
and make  sure the paths containing your .jar libraries are listed. For example, 
my ```classpath.txt``` file includes the path where the java compiled sofware of JDemetra+ is included:
```L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB\demetra-tstoolkit-2.2.2.jar```
If you don't want to modify the ```classpath.txt``` file because you are using sofware that relies on older Java versions, then            add the line:
```javaclasspath('L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB\')``` 
at the beginning of the ```mjdemetra``` function so that the desired version of Java is used only within within the function.

## Example

The ```mjdemetra``` function can be used in many different ways. By default it plots the seasonally adjusted data (with confidence intervals only when the TramoSeats method is used) and highlights outliers. The seasonally adjusted series without removing calendar effects is also plotted.

```
        [sa, rslts]= mjdemetra2(data,'horizon',20,'Method','TramoSeats','CalendarOption','RSAfull')
        [sa, rslts]= mjdemetra2(data2,            'Method','X13'      );
        [sa, rslts]= mjdemetra2(data,'horizon',20,'Method','TramoSeats','CalendarOption','RSA5')
        [sa, rslts]= mjdemetra2(data,'horizon',20,'Method','X13'       ,'CalendarOption','RSA5c')
        [sa, rslts]= mjdemetra2(data)
        [sa, rslts]= mjdemetra2(data,                                  ,'CalendarOption','RSA0')
        [sa, rslts]= mjdemetra2(data,                                                          , 'grafico',false)
```
