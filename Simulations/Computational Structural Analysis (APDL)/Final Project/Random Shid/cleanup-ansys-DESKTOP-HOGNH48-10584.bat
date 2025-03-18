@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 8836)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 6360)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 6232)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 10584)

del /F cleanup-ansys-DESKTOP-HOGNH48-10584.bat
