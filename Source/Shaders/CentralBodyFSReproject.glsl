uniform sampler2D u_texture;

uniform vec2 u_south;
uniform vec2 u_deltaWLat;
uniform vec2 u_minMLat;
uniform vec2 u_maxMLat;

uniform float u_height;

varying vec2 v_textureCoordinates;

void main()
{
    vec2 texCoords = v_textureCoordinates;
    if (texCoords.y > 0.0 && texCoords.y < 1.0)
    {
        vec2 currentWLat = agi_df64Add(u_south, agi_df64Multiply(u_deltaWLat, texCoords.y * u_height));
        
        vec2 sinTheta = agi_df64Sin2(currentWLat);
        vec2 mLat = agi_df64Multiply(0.5, agi_df64Log(agi_df64Divide(agi_df64Add(1.0, sinTheta), agi_df64Subtract(1.0, sinTheta))));
        vec2 mRow = agi_df64Divide(agi_df64Subtract(mLat, u_minMLat), agi_df64Subtract(u_maxMLat, u_minMLat));
        
        texCoords.y = (mRow.x + mRow.y);
    }
    
    gl_FragColor = texture2D(u_texture, texCoords);
}
