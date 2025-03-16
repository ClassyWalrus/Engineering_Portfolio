@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 4624)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 2676)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 13832)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 13156)

del /F cleanup-ansys-DESKTOP-HOGNH48-13156.bat
