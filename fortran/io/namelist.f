namelist /meeting/joe,jake,jane
        open(1,file='namefile')
        read(1,meeting)
        print *, 'joe = ',joe,' jake = ',jake,' jane = ',jane
        end

