uniform sampler2D u_texture;

uniform float u_north;
uniform float u_south;
uniform float u_width;
uniform float u_height;

varying vec2 v_textureCoordinates;

void main()
{
    vec2 texCoords = v_textureCoordinates;
    if (texCoords.y > 0.0 && texCoords.y < 1.0)
    {
        float sinTheta = sin(u_south);
        float minMLat = 0.5 * log((1.0 + sinTheta) / (1.0 - sinTheta));
        sinTheta = sin(u_north);
        float maxMLat = 0.5 * log((1.0 + sinTheta) / (1.0 - sinTheta));
        float invMLatDim = 1.0 / (maxMLat - minMLat);
        
        float heightMinusOne = u_height - 1.0;
        float deltaWLat = (u_north - u_south) / u_height;
        float currentWLat = u_south + deltaWLat * (0.5 + texCoords.y * u_height);
        
        sinTheta = sin(currentWLat);
        float mLat = 0.5 * log((1.0 + sinTheta) / (1.0 - sinTheta));
        float mRow = floor(heightMinusOne * (mLat - minMLat) * invMLatDim);
        
        texCoords.y = mRow / u_height;
    }
    
    gl_FragColor = texture2D(u_texture, texCoords);
}