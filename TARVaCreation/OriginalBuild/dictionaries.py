class Dictionaries:
    @staticmethod

    def _one_key(d, k1, kind):
        d[k1] = kind
        return(d)

    def _one_key_update(d,k1,kind,u):
        d[k1] = kind
        if type(kind) == list:
            d[k1].append(u)
        else:
            d[k1] = u
        return(d)

    def _one_just_update(d,k1,u):
        if type(d[k1]) == list:
            d[k1].append(u)
        else:
            d[k1] = u
        return(d)

    def _two_key(d,k1,k2,kind):
        d[k1][k2] = kind
        return(d)
    
    def _two_key_update(d, k1,k2,kind,u):
        d[k1][k2] = kind
        if type(kind) == list:
            d[k1][k2].append(u)
        else:
            d[k1][k2] = u
        return(d)

    def _two_just_update(d,k1,k2,u):
        if type(d[k1][k2]) == list:
            d[k1][k2].append(u)
        else:
            d[k1][k2] = u
        return(d)

    def _three_key(d,k1,k2,k3,kind):
        d[k1][k2][k3] = kind
        return(d)

    def _three_key_update(d,k1,k2,k3,kind,u):
        d[k1][k2][k3] = kind
        if type(d[k1][k2][k3]) == list:
            d[k1][k2][k3].append(u)
        else:
            d[k1][k2][k3] = u
        return(d)

    def _three_just_update(d,k1,k2,k3,u):
        if type(d[k1][k2][k3]) == list:
            d[k1][k2][k3].append(u)
        else:
            d[k1][k2][k3] = u
        return(d)

    def _four_key(d,k1,k2,k3,k4,kind):
        d[k1][k2][k3][k4] = kind
        return(d)

    def _four_key_update(d,k1,k2,k3,k4,kind,u):
        d[k1][k2][k3][k4] = kind
        if type(d[k1][k2][k3][k4]) == list:
            d[k1][k2][k3][k4].append(u)
        else:
            d[k1][k2][k3][k4] = u
        return(d)

    def _four_just_update(d,k1,k2,k3,k4,u):
        if type(d[k1][k2][k3][k4]) == list:
            d[k1][k2][k3][k4].append(u)
        else:
            d[k1][k2][k3][k4] = u
        return(d)

    def _five_key(d,k1,k2,k3,k4,k5,kind):
        d[k1][k2][k3][k4][k5] = kind
        return(d)

    def _five_key_update(d,k1,k2,k3,k4,k5,kind,u):
        d[k1][k2][k3][k4][k5] = kind
        if type(d[k1][k2][k3][k4][k5]) == list:
            d[k1][k2][k3][k4][k5].append(u)
        else:
            d[k1][k2][k3][k4][k5] = u
        return(d)

    def _five_just_update(d,k1,k2,k3,k4,k5,u):
        if type(d[k1][k2][k3][k4][k5]) == list:
            d[k1][k2][k3][k4][k5].append(u)
        else:
            d[k1][k2][k3][k4][k5] = u
        return(d)
    
    def _six_new_kv_pairs(d,k1,k2,k3,k4,k5,k6,v1,v2,v3,v4,v5,v6):
        d[k1] = v1
        d[k2] = v2
        d[k3] = v3
        d[k4] = v4
        d[k5] = v5
        d[k6] = v6
        return(d)

    def _loop_list_kvs(d,k_list,v_list):
        for k in range(0, len(k_list)):
            key = k_list[k]
            valu = v_list[k]
            d[key] = valu
        return(d)

    

