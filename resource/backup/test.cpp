        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          if (!inv_bhead[i]){
            // Non-basic
            if (isGreater(alpha_r(i), 0) && isEqual(x(i), l(i))){
              local_init_size ++;
              local_init.push_back(i);
            }
          }
        }
        #pragma omp for nowait
        for (int i = 0; i < n; i ++){
          if (!inv_bhead[i]){
            // Non-basic
            if (isLess(alpha_r(i), 0) && isEqual(x(i), u(i))){
              local_init_size ++;
              local_init.push_back(i);
            }
          }
        }
        #pragma omp single nowait
        {
          for (int i = n; i < n+m; i ++){
            if (!inv_bhead[i]){
              // Non-basic
              if ((isGreater(alpha_r(i), 0) && isEqual(x(i), bl(i-n))) || (isLess(alpha_r(i), 0) && isEqual(x(i), bu(i-n)))){
                local_init_size ++;
                local_init.push_back(i);
              }
            }
          }
        }