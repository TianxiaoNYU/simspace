import numpy as np

def create_ellipse(
        array_shape: tuple, 
        center: tuple, 
        radius_x: int, 
        radius_y: int, 
        angle: float) -> np.ndarray:
    """
    Create a 2D numpy array with a rotated ellipse shape.

    Args:
        array_shape (tuple): Shape of the 2D array.
        center (tuple): Center of the ellipse as a tuple (x, y).
        radius_x (int): Radius of the ellipse along the x-axis.
        radius_y (int): Radius of the ellipse along the y-axis.
        angle (float): Angle of rotation in degrees.

    Returns:
        numpy.ndarray: 2D numpy array with the rotated ellipse shape, where pixels inside the ellipse are set to 1 and others to 0.

    Examples:
        >>> import numpy as np
        >>> arr = create_ellipse((10, 10), (5, 5), 3, 2, 45)
        >>> print(arr.astype(int))
        [[0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 1 0 0 0 0]
        [0 0 0 1 1 1 1 0 0 0]
        [0 0 1 1 1 1 1 1 0 0]
        [0 0 1 1 1 1 1 1 0 0]
        [0 0 0 1 1 1 1 0 0 0]
        [0 0 0 0 1 1 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]]   
    """
    # Create an empty numpy array with the given shape
    array = np.zeros(array_shape)
    
    # Get the x and y coordinates of the center
    center_x, center_y = center
    
    # Convert the angle from degrees to radians
    angle_rad = np.radians(angle)
    
    # Iterate over each pixel in the array
    for x in range(array_shape[0]):
        for y in range(array_shape[1]):
            # Calculate the coordinates of the current pixel relative to the center of the ellipse
            x_rel = x - center_x
            y_rel = y - center_y
            
            # Apply the rotation transformation
            x_rot = x_rel * np.cos(angle_rad) - y_rel * np.sin(angle_rad)
            y_rot = x_rel * np.sin(angle_rad) + y_rel * np.cos(angle_rad)
            
            # Calculate the distance from the current pixel to the center of the rotated ellipse
            distance_squared = ((x_rot / radius_x) ** 2) + ((y_rot / radius_y) ** 2)
            
            # If the distance is less than or equal to 1, set the pixel value to 1
            if distance_squared <= 1:
                array[x, y] = 1
    
    return array

def create_vessel(array_shape, center, length, width):
    """
    Create a 2D numpy array with a linear tube shape.

    Args:
        array_shape (tuple): Shape of the 2D array.
        center (tuple): Center of the tube as a tuple (x, y).
        length (int): Length of the tube.
        width (int): Width of the tube.

    Returns:
        numpy.ndarray: 2D numpy array with the tube shape

    Examples:
        >>> import numpy as np
        >>> arr = create_vessel((10, 10), (5, 5), 6, 2)
        >>> print(arr.astype(int))
        [[0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 1 1 1 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]]]
    """
    # Create an empty numpy array with the given shape
    array = np.zeros(array_shape)
        
    # Get the x and y coordinates of the center
    center_x, center_y = center
        
    # Calculate the starting and ending coordinates of the tube
    start_x = center_x - length // 2
    end_x = center_x + length // 2
    start_y = center_y - width // 2
    end_y = center_y + width // 2
        
    # Iterate over each pixel in the array
    for x in range(array_shape[0]):
        for y in range(array_shape[1]):
            # If the pixel is within the tube's boundaries, set the pixel value to 1
            if start_x <= x <= end_x and start_y <= y <= end_y:
                array[x, y] = 1
        
    return array

def create_ring(array_shape, center, inner_radius, outer_radius):
    """
    Create a 2D numpy array with a ring shape.

    Args:
        array_shape (tuple): Shape of the 2D array.
        center (tuple): Center of the ring as a tuple (x, y).
        inner_radius (int): Inner radius of the ring.
        outer_radius (int): Outer radius of the ring.

    Returns:
        numpy.ndarray: 2D numpy array with the ring shape

    Examples:
        >>> import numpy as np
        >>> arr = create_ring((10, 10), (5, 5), 2, 4)
        >>> print(arr.astype(int))
        [[0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 1 0 0 0 0]
        [0 0 0 1 1 1 1 1 0 0]
        [0 0 1 1 1 1 1 1 1 0]
        [0 0 1 1 0 0 0 1 1 0]
        [0 1 1 1 0 0 0 1 1 1]
        [0 0 1 1 0 0 0 1 1 0]
        [0 0 1 1 1 1 1 1 1 0]
        [0 0 0 1 1 1 1 1 0 0]
        [0 0 0 0 0 1 0 0 0 0]]
    """
    # Create an empty numpy array with the given shape
    array = np.zeros(array_shape)
    
    # Get the x and y coordinates of the center
    center_x, center_y = center
    
    # Calculate the squared radii
    inner_radius_squared = inner_radius ** 2
    outer_radius_squared = outer_radius ** 2
    
    # Iterate over each pixel in the array
    for x in range(array_shape[0]):
        for y in range(array_shape[1]):
            # Calculate the distance from the current pixel to the center of the ring
            distance_squared = (x - center_x) ** 2 + (y - center_y) ** 2
            
            # If the distance is within the ring boundaries, set the pixel value to 1
            if inner_radius_squared <= distance_squared <= outer_radius_squared:
                array[x, y] = 1
    
    return array

def create_rectangle(array_shape, center, width, height):
    """
    Create a 2D numpy array with a rectangle shape.

    Args:
        array_shape (tuple): Shape of the 2D array.
        center (tuple): Center of the rectangle as a tuple (x, y).
        width (int): Width of the rectangle.
        height (int): Height of the rectangle.

    Returns:
        numpy.ndarray: 2D numpy array with the rectangle shape

    Examples:
        >>> import numpy as np
        >>> arr = create_rectangle((10, 10), (5, 5), 2, 4)
        >>> print(arr.astype(int))
        [[0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 1 1 1 1 1 0 0]
        [0 0 0 1 1 1 1 1 0 0]
        [0 0 0 1 1 1 1 1 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]
        [0 0 0 0 0 0 0 0 0 0]]
    """
    # Create an empty numpy array with the given shape
    array = np.zeros(array_shape)
    
    # Get the x and y coordinates of the center
    center_x, center_y = center
    
    # Calculate the starting and ending coordinates of the rectangle
    start_x = center_x - width // 2
    end_x = center_x + width // 2
    start_y = center_y - height // 2
    end_y = center_y + height // 2
    
    # Iterate over each pixel in the array
    for x in range(array_shape[0]):
        for y in range(array_shape[1]):
            # If the pixel is within the rectangle's boundaries, set the pixel value to 1
            if start_x <= x <= end_x and start_y <= y <= end_y:
                array[x, y] = 1
    
    return array