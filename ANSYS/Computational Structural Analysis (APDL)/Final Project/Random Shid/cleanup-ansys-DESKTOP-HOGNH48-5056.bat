@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 20380)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14584)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 20148)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 5056)

del /F cleanup-ansys-DESKTOP-HOGNH48-5056.bat
