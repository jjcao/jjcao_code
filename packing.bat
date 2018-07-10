@Echo Off

del *.ncb /s
del *.sdf /s
del *.VC.db /s
del *.asv /s
del UpgradeLog.XML /s

@Echo Find debug
@For /r . %%a In (.) Do @If Exist "%%a\debug" @Echo "%%a\debug"

@Echo Find debug Dir....OK
@For /r . %%a In (.) Do @If Exist "%%a\debug" rd /s /q "%%a\debug"
@Echo Clear debug Dir Mission Completed


@Echo Find release
@For /r . %%a In (.) Do @If Exist "%%a\release" @Echo "%%a\release"

@Echo Find release Dir....OK
@For /r . %%a In (.) Do @If Exist "%%a\release" rd /s /q "%%a\release"
@Echo Clear rease Dir Mission Completed


@Echo Find ipch
@For /r . %%a In (.) Do @If Exist "%%a\ipch" @Echo "%%a\ipch"
@Echo Find ipch Dir....OK
@For /r . %%a In (.) Do @If Exist "%%a\ipch" rd /s /q "%%a\ipch"


@Echo Find ipch
@For /r . %%a In (.) Do @If Exist "%%a\_UpgradeReport_Files" @Echo "%%a\_UpgradeReport_Files"
@For /r . %%a In (.) Do @If Exist "%%a\_UpgradeReport_Files" rd /s /q "%%a\_UpgradeReport_Files"


@echo on
