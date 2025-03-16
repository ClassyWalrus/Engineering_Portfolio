@echo off
set LOCALHOST=%COMPUTERNAME%
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 19792)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 10376)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 5148)
if /i "%LOCALHOST%"=="DESKTOP-HOGNH48" (taskkill /f /pid 14616)

del /F cleanup-ansys-DESKTOP-HOGNH48-14616.bat
