function response = argoGeminiGenerateContent(prompt, nvp)
    arguments
        prompt      (1,1) {mustBeTextScalar,mustBeNonempty}
        nvp.image   (1,1) {mustBeTextScalar} = "";
    end

    if strlength(nvp.image) == 0
        model = "argo:gemini-2.5-flash";
        query = struct();
        query.model = model;
        query.prompt = prompt;
        query.user = 'mmenickelly';
    else
        % THIS WILL NOT WORK AS IS. 
        model = "gemini-pro-vision";
        if contains(nvp.image,["http://","https://"])
            imdata = imread(nvp.image);
            imwrite(imdata,"imdata.png")
            img = "imdata.png";
        else
            img = nvp.image;
        end
        fid = fopen(img);
        im = fread(fid,'*uint8');
        fclose(fid);
        b64 = matlab.net.base64encode(im);
        [~,~,ext] = fileparts(img);
        MIMEType = "image/" + erase(ext,".");
        query = struct("contents",[]);
        query.contents = {struct("parts",[])};
        query.contents{1}.parts = {struct("text",prompt),struct("inline_data",[])};
        query.contents{1}.parts{2}.inline_data = struct("mime_type",MIMEType,"data",[]);
        query.contents{1}.parts{2}.inline_data.data = b64;
        if isfile("imdata.png")
            delete("imdata.png")
        end
    end
    endpoint = "http://localhost:63332/v1/chat";
    
    import matlab.net.*
    import matlab.net.http.*  
    headers = HeaderField('Content-Type', 'application/json');
    request = RequestMessage('post', headers, query);
    response = send(request, URI(endpoint));
end