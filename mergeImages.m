image_path = 'C:\Users\Peihong\Desktop\Code\Utilities\LaplacianDeformation\clips\1\poses\';

for i=1:18
    It = imread([image_path, 'target ', num2str(i), '.png']);
    Id = imread([image_path, 'deformed ', num2str(i), '.png']);
    Ie = imread([image_path, 'error ', num2str(i), '.png']);
    Io = imread([image_path, 'overlay ', num2str(i), '.png']);
    
    I = [It, Id; Ie, Io];
    imwrite(I, [image_path, 'combined ', num2str(i), '.png']);
end