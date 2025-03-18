@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 15364)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 10392)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14128)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 19152)

del /F cleanup-ansys-DESKTOP-HOGNH48-19152.bat
