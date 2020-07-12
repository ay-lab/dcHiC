setlocal
::Get the current batch file's short path
for %%x in (%0) do set BatchPath=%%~dpsx
for %%x in (%BatchPath%) do set BatchPath=%%~dpsx

if exist %BatchPath%\jdk-11 (
  echo "Using bundled JDK."
  set JAVA_HOME=%BatchPath%\jdk-11
  set JAVA_CMD=%BatchPath%\jdk-11\bin\javaw
) else (
  echo "Bundled JDK not found.  Using system JDK."
  set JAVA_CMD=java
)

start %JAVA_CMD% -showversion --module-path=%BatchPath%\modules -Xmx4g -Djava.awt.headless=true @%BatchPath%\gsea.args --patch-module=jide.common=%BatchPath%\lib\jide-components-3.7.4.jar;%BatchPath%\lib\jide-dock-3.7.4.jar;%BatchPath%\lib\jide-grids-3.7.4.jar --module=org.gsea_msigdb.gsea/xapps.gsea.CLI  %*
