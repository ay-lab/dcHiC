#!/bin/sh

#This script is intended for launch on *nix machines

#-Xmx4g indicates 4 gb of memory, adjust number up or down as needed
prefix=`dirname $(readlink $0 || echo $0)`

# Check whether or not to use the bundled JDK
if [ -d "${prefix}/jdk-11" ]; then
    echo echo "Using bundled JDK."
    JAVA_HOME="${prefix}/jdk-11"
    PATH=$JAVA_HOME/bin:$PATH
else
    echo "Using system JDK."
fi

if [ -e "${prefix}/modules/disable-prefs.jar" ]; then
    # Running in a context with Preferences disabled (probably as a GP Module)
    PREFS_PROP=-Djava.util.prefs.PreferencesFactory=com.allaboutbalance.articles.disableprefs.DisabledPreferencesFactory
else
    PREFS_PROP=
fi;

exec java --module-path="${prefix}/modules" -Xmx4g \
    -Djava.awt.headless=true $PREFS_PROP \
    @"${prefix}/gsea.args" \
    --patch-module="jide.common=${prefix}/lib/jide-components-3.7.4.jar:${prefix}/lib/jide-dock-3.7.4.jar:${prefix}/lib/jide-grids-3.7.4.jar" \
    --module=org.gsea_msigdb.gsea/xapps.gsea.CLI "$@"
