@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 3872)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 27428)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 3580)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 24208)

del /F cleanup-ansys-DESKTOP-HOGNH48-24208.bat
