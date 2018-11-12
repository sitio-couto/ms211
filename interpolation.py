x = [1.7,2,2.5,3]
y = [-4.2,2,-1.4,0.8]

def print_diff_table(orders):
    for i in range(len(orders)):
        print("Order", i)
        for x in orders[i]:
            print("%.2f" % x)
        print("\n")

def diff_table(x, y):
    orders = [y.copy()]

    for ord in range(1,len(x)):
        new_ord = []
        for i in range(len(x)-ord):
            last_ord = orders[-1]
            new_ord.append((last_ord[i+1]-last_ord[i])/(x[i+ord]-x[i]))
        orders.append(new_ord.copy())

    print_diff_table(orders)
    return orders


orders = diff_table(x, y)
