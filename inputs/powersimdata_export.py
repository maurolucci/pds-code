from powersimdata import Grid

if __name__ == '__main__':
    interconnects = ['Eastern', 'Western', 'Texas']
    for con in interconnects:
        grid = Grid(con)
        print(len(grid.bus))
        print(len(grid.branch))
        print(dict(grid.plant))