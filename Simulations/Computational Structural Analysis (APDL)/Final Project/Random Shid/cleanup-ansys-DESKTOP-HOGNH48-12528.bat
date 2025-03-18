@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 18120)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 6728)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 10712)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 12528)

del /F cleanup-ansys-DESKTOP-HOGNH48-12528.bat
