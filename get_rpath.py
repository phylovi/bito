import site
print('-Wl,-rpath,.,-rpath,' + site.getusersitepackages() + '/bito,-rpath,' + ',-rpath,'.join(site.getsitepackages()) + '/bito')
